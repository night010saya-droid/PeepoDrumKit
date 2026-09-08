#include "audio_tempo_analysis.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace Audio
{
	namespace
	{
		struct Onset
		{
			f64 Frame = 0.0;
			f32 Strength = 0.0f;
		};

		constexpr f64 MinBPM = 90.0;
		constexpr f64 MaxBPM = 205.0;
		constexpr f64 CoarseBPMStep = 0.5;
		constexpr f64 FineBPMStep = 0.001;
		constexpr i64 AnalysisWindow = 1024;
		constexpr i64 AnalysisHop = 512;
		constexpr i64 MinimumOnsetDistanceMS = 80;
		constexpr i64 GapWindowSamples = 2048;

		static std::vector<Onset> FindOnsets(const PCMSampleBuffer& buffer, i64 firstFrame, i64 frameCount)
		{
			const i64 lastFrame = firstFrame + frameCount;
			const i64 energyBegin = std::max<i64>(0, firstFrame - AnalysisHop - AnalysisWindow);
			const i64 energyFrameCount = lastFrame - energyBegin;
			std::vector<f64> energyPrefix(static_cast<size_t>(energyFrameCount + 1), 0.0);
			for (i64 frame = 0; frame < energyFrameCount; ++frame)
			{
				f64 mixed = 0.0;
				const i64 sourceFrame = energyBegin + frame;
				for (u32 c = 0; c < buffer.ChannelCount; ++c)
					mixed += ConvertSampleI16ToF32(buffer.InterleavedSamples[sourceFrame * buffer.ChannelCount + c]);
				mixed /= static_cast<f64>(buffer.ChannelCount);
				energyPrefix[static_cast<size_t>(frame + 1)] = energyPrefix[static_cast<size_t>(frame)] + mixed * mixed;
			}

			const auto frameRMS = [&](i64 frame, i64 count)
			{
				const i64 begin = Clamp(frame, energyBegin, lastFrame);
				const i64 end = Clamp(frame + count, begin, lastFrame);
				if (end <= begin)
					return 0.0f;
				const f64 sum = energyPrefix[static_cast<size_t>(end - energyBegin)] - energyPrefix[static_cast<size_t>(begin - energyBegin)];
				return static_cast<f32>(std::sqrt(sum / static_cast<f64>(end - begin)));
			};

			std::vector<f32> flux;
			for (i64 frame = firstFrame; frame < lastFrame; frame += AnalysisHop)
			{
				const f32 current = frameRMS(frame, AnalysisWindow);
				const f32 previous = frameRMS(frame - AnalysisHop, AnalysisWindow);
				flux.push_back(std::max(0.0f, current - previous));
			}

			if (flux.size() < 3)
				return {};

			const f32 mean = std::accumulate(flux.begin(), flux.end(), 0.0f) / static_cast<f32>(flux.size());
			f32 variance = 0.0f;
			for (f32 value : flux)
				variance += (value - mean) * (value - mean);
			const f32 threshold = mean + std::sqrt(variance / static_cast<f32>(flux.size())) * 0.45f;
			const i64 minimumDistance = std::max<i64>(1, static_cast<i64>(buffer.SampleRate) * MinimumOnsetDistanceMS / 1000 / AnalysisHop);

			std::vector<Onset> onsets;
			for (size_t i = 1; i + 1 < flux.size(); ++i)
			{
				if (flux[i] < threshold || flux[i] < flux[i - 1] || flux[i] < flux[i + 1])
					continue;

				f64 peakOffset = 0.0;
				const f64 curvature = static_cast<f64>(flux[i - 1]) - 2.0 * static_cast<f64>(flux[i]) + static_cast<f64>(flux[i + 1]);
				if (std::abs(curvature) > 0.000001)
					peakOffset = Clamp(0.5 * (static_cast<f64>(flux[i - 1]) - static_cast<f64>(flux[i + 1])) / curvature, -0.5, 0.5);
				const Onset onset { static_cast<f64>(firstFrame + static_cast<i64>(i) * AnalysisHop) + peakOffset * AnalysisHop, flux[i] };
				if (!onsets.empty() && onset.Frame - onsets.back().Frame < minimumDistance * AnalysisHop)
				{
					if (onset.Strength > onsets.back().Strength)
						onsets.back() = onset;
				}
				else
					onsets.push_back(onset);
			}
			return onsets;
		}

		static f64 CircularWindowSupport(const std::vector<f64>& doubledPrefix, i64 histogramSize, i64 center, i64 windowSize)
		{
			if (histogramSize <= 0)
				return 0.0;
			const i64 halfWindow = std::min<i64>(windowSize / 2, (histogramSize - 1) / 2);
			const i64 normalizedCenter = (center % histogramSize + histogramSize) % histogramSize;
			i64 begin = normalizedCenter - halfWindow;
			if (begin < 0)
				begin += histogramSize;
			const i64 end = begin + halfWindow * 2 + 1;
			if (begin < 0 || end < begin || end >= static_cast<i64>(doubledPrefix.size()))
				return 0.0;
			return doubledPrefix[static_cast<size_t>(end)] - doubledPrefix[static_cast<size_t>(begin)];
		}

		static TempoAnalysisCandidate EvaluateInterval(const std::vector<Onset>& onsets, u32 sampleRate, i64 interval, i32 downsample)
		{
			const i64 reducedInterval = std::max<i64>(1, interval >> downsample);
			const i64 sampleStep = i64(1) << downsample;
			std::vector<f64> histogram(static_cast<size_t>(reducedInterval), 0.0);
			f64 totalStrength = 0.0;
			for (const Onset& onset : onsets)
			{
				const i64 phase = static_cast<i64>(std::fmod(onset.Frame, static_cast<f64>(interval)));
				const i64 position = Clamp<i64>(phase / sampleStep, 0, reducedInterval - 1);
				histogram[static_cast<size_t>(position)] += onset.Strength;
				totalStrength += onset.Strength;
			}
			std::vector<f64> doubledPrefix(static_cast<size_t>(reducedInterval * 2 + 1), 0.0);
			for (i64 i = 0; i < reducedInterval * 2; ++i)
				doubledPrefix[static_cast<size_t>(i + 1)] = doubledPrefix[static_cast<size_t>(i)] + histogram[static_cast<size_t>(i % reducedInterval)];

			const i64 windowSize = std::min<i64>(reducedInterval, std::max<i64>(8, GapWindowSamples >> downsample));
			f64 bestScore = 0.0;
			i64 bestPosition = 0;
			for (const Onset& onset : onsets)
			{
				const i64 phase = static_cast<i64>(std::fmod(onset.Frame, static_cast<f64>(interval)));
				const i64 position = Clamp<i64>(phase / sampleStep, 0, reducedInterval - 1);
				const i64 offbeatPosition = (position + reducedInterval / 2) % reducedInterval;
				const f64 score = CircularWindowSupport(doubledPrefix, reducedInterval, position, windowSize) + CircularWindowSupport(doubledPrefix, reducedInterval, offbeatPosition, windowSize) * 0.5;
				if (score > bestScore)
				{
					bestScore = score;
					bestPosition = position;
				}
			}

			TempoAnalysisCandidate result;
			result.BPM = static_cast<f32>(static_cast<f64>(sampleRate) * 60.0 / static_cast<f64>(interval));
			result.Offset = Time::FromSec(static_cast<f64>(bestPosition * sampleStep) / static_cast<f64>(sampleRate));
			result.Confidence = (totalStrength > 0.0) ? static_cast<f32>(bestScore / totalStrength) : 0.0f;
			return result;
		}

		static std::vector<f32> ComputeWaveformSlopes(const PCMSampleBuffer& buffer)
		{
			std::vector<f32> slopes(static_cast<size_t>(buffer.FrameCount), 0.0f);
			const i64 halfWindow = std::max<i64>(1, buffer.SampleRate / 20);
			if (buffer.FrameCount < halfWindow * 2)
				return slopes;

			const auto absoluteSample = [&](i64 frame)
			{
				f64 mixed = 0.0;
				for (u32 c = 0; c < buffer.ChannelCount; ++c)
					mixed += std::abs(static_cast<f64>(ConvertSampleI16ToF32(buffer.InterleavedSamples[frame * buffer.ChannelCount + c])));
				return mixed / static_cast<f64>(buffer.ChannelCount);
			};

			f64 left = 0.0;
			f64 right = 0.0;
			for (i64 i = 0; i < halfWindow; ++i)
			{
				left += absoluteSample(i);
				right += absoluteSample(i + halfWindow);
			}
			for (i64 frame = halfWindow; frame < buffer.FrameCount - halfWindow; ++frame)
			{
				slopes[static_cast<size_t>(frame)] = static_cast<f32>(std::max(0.0, (right - left) / static_cast<f64>(halfWindow)));
				const f64 current = absoluteSample(frame);
				left += current - absoluteSample(frame - halfWindow);
				right += absoluteSample(frame + halfWindow) - current;
			}
			return slopes;
		}

		static f64 SampleSlope(const std::vector<f32>& slopes, f64 frame)
		{
			if (slopes.empty() || frame < 0.0 || frame >= static_cast<f64>(slopes.size() - 1))
				return 0.0;
			const i64 index = static_cast<i64>(frame);
			const f64 t = frame - static_cast<f64>(index);
			return slopes[static_cast<size_t>(index)] * (1.0 - t) + slopes[static_cast<size_t>(index + 1)] * t;
		}

		static void RefineOffsetWithWaveform(const std::vector<f32>& slopes, u32 sampleRate, TempoAnalysisCandidate& candidate)
		{
			const f64 secondsPerBeat = 60.0 / candidate.BPM;
			const f64 interval = secondsPerBeat * sampleRate;
			f64 offset = candidate.Offset.ToSec() * sampleRate;
			const f64 offbeat = std::fmod(offset + interval * 0.5, interval);
			f64 scoreA = 0.0;
			f64 scoreB = 0.0;
			for (f64 positionA = offset, positionB = offbeat; positionA < slopes.size() && positionB < slopes.size(); positionA += interval, positionB += interval)
			{
				scoreA += SampleSlope(slopes, positionA);
				scoreB += SampleSlope(slopes, positionB);
			}
			if (scoreB > scoreA)
				offset = offbeat;
			candidate.Offset = Time::FromMS(std::round(Time::FromSec(offset / sampleRate).ToMS()));
		}
	}

	TempoAnalysisResult AnalyzeTempo(const PCMSampleBuffer& buffer, Time start, Time duration)
	{
		TempoAnalysisResult result;
		if (buffer.ChannelCount == 0 || buffer.SampleRate == 0 || buffer.FrameCount <= AnalysisWindow)
			return result;

		const i64 firstFrame = Clamp<i64>(TimeToFrames(start, buffer.SampleRate), 0, buffer.FrameCount - 1);
		const i64 availableFrames = buffer.FrameCount - firstFrame;
		const i64 requestedFrames = (duration > Time::Zero()) ? TimeToFrames(duration, buffer.SampleRate) : availableFrames;
		const i64 frameCount = Clamp<i64>(requestedFrames, AnalysisWindow, availableFrames);
		result.AnalyzedDuration = FramesToTime(frameCount, buffer.SampleRate);

		const std::vector<Onset> onsets = FindOnsets(buffer, firstFrame, frameCount);
		if (onsets.size() < 2)
			return result;

		const i64 minInterval = static_cast<i64>(std::ceil(static_cast<f64>(buffer.SampleRate) * 60.0 / MaxBPM));
		const i64 maxInterval = static_cast<i64>(std::floor(static_cast<f64>(buffer.SampleRate) * 60.0 / MinBPM));
		std::vector<TempoAnalysisCandidate> candidates;
		for (i64 interval = minInterval; interval <= maxInterval; interval += 10)
			candidates.push_back(EvaluateInterval(onsets, buffer.SampleRate, interval, 3));
		std::stable_sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.Confidence > b.Confidence; });

		// Refine several promising sample intervals at full resolution. This is
		// the same interval-first search strategy used by ArrowVortex.
		const size_t coarseCandidateCount = std::min<size_t>(candidates.size(), 8);
		std::vector<TempoAnalysisCandidate> refinedCandidates;
		for (size_t i = 0; i < coarseCandidateCount; ++i)
		{
			const i64 centerInterval = static_cast<i64>(std::llround(static_cast<f64>(buffer.SampleRate) * 60.0 / candidates[i].BPM));
			const i64 begin = std::max(minInterval, centerInterval - 10);
			const i64 end = std::min(maxInterval, centerInterval + 10);
			for (i64 interval = begin; interval <= end; ++interval)
				refinedCandidates.push_back(EvaluateInterval(onsets, buffer.SampleRate, interval, 0));
		}
		std::stable_sort(refinedCandidates.begin(), refinedCandidates.end(), [](const auto& a, const auto& b) { return a.Confidence > b.Confidence; });

		for (const TempoAnalysisCandidate& candidate : refinedCandidates)
		{
			if (result.CandidateCount >= result.Candidates.size())
				break;
			b8 duplicate = false;
			for (size_t i = 0; i < result.CandidateCount; ++i)
				duplicate |= Absolute(candidate.BPM - result.Candidates[i].BPM) < 4.0f;
			if (!duplicate)
				result.Candidates[result.CandidateCount++] = candidate;
		}

		if (result.IsValid())
		{
			const std::vector<f32> slopes = ComputeWaveformSlopes(buffer);
			for (size_t i = 0; i < result.CandidateCount; ++i)
				RefineOffsetWithWaveform(slopes, buffer.SampleRate, result.Candidates[i]);
		}
		return result;
	}
}
