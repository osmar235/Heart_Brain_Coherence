@echo off
setlocal EnableExtensions
rem Step 09: group statistics exactly as frozen in docs\ANALYSIS_FREEZE_v1.md. Seconds to run.
rem Run it once, after step 08 has finished for the whole cohort.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.cohort_stats %*
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
:fail
popd
exit /b 1
