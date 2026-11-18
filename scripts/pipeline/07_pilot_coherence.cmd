@echo off
setlocal EnableExtensions
rem Step 07: first coherence run on the pilot. FEASIBILITY ONLY (docs\PILOT_COHERENCE_PLAN.md).
rem Needs step 03b and step 06 outputs. Reads no raw recording. A few minutes.
rem No span against span+1 comparison and no p values are produced.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.pilot_coherence %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo STEP 07 DONE. See results\pilot_coherence
if not "%RC%"=="0" echo STEP 07 FAILED, exit code %RC%. Send back the newest logs\pilot_coherence_*.log
popd
exit /b %RC%
:fail
popd
exit /b 1
