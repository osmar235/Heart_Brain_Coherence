@echo off
setlocal EnableExtensions
rem Step 03b: applies the proposed beat correction policy (D-007) to the nine pilot recordings.
rem Reads each recording once (a few minutes, about 4 GB of free memory). Changes no raw data.
rem Writes results\beat_correction and derivatives\beats.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.beat_correction %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo STEP 03b DONE. See results\beat_correction
if not "%RC%"=="0" echo STEP 03b NOT CLEAN, exit code %RC%. Send back the newest logs\beat_correction_*.log
popd
exit /b %RC%
:fail
popd
exit /b 1
