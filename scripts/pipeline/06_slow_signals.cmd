@echo off
setlocal EnableExtensions
rem Step 06: slow 4 Hz heart and EEG envelope series for the pilot, plus trial locked averages.
rem Needs step 03b outputs (derivatives\beats). Reads each recording once; about 10 minutes in total
rem and about 5 GB of free memory. EEG handling here is provisional (see src\hb\slow_signals.py).
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.slow_signals %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo STEP 06 DONE. See results\slow_signals
if not "%RC%"=="0" echo STEP 06 NOT CLEAN, exit code %RC%. Send back the newest logs\slow_signals_*.log
popd
exit /b %RC%
:fail
popd
exit /b 1
