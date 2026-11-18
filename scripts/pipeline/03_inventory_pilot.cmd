@echo off
setlocal EnableExtensions
rem Step 03: inventory and ECG quality check of the pilot recordings.
rem Reads the nine pilot files one at a time. Needs about 4 GB of free memory and a few minutes.
rem Writes results\pilot_qc (figures and tables) and derivatives\rpeaks. Changes no raw data.
rem Optional: --subjects sub-124 --tasks video rest
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -c "import numpy, scipy, matplotlib" 
if errorlevel 1 goto :missing
"%PY%" -W ignore -m hb.pilot_qc %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo STEP 03 DONE. See results\pilot_qc
if not "%RC%"=="0" echo STEP 03 NOT CLEAN, exit code %RC%. Send back the newest logs\pilot_qc_*.log
popd
exit /b %RC%
:missing
echo A required package is missing. Install with:  "%PY%" -m pip install -r requirements.txt
:fail
popd
exit /b 1
