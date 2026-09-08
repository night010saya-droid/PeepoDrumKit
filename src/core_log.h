#pragma once

#include "core_types.h"

namespace Log
{
	void Initialize();
	void Shutdown();
	void Write(cstr format, ...);
	void WriteHRESULT(cstr operation, long result);
	void InstallCrashHandler();
}
