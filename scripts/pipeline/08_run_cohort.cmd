@echo off
setlocal EnableExtensions
rem Step 08: whole cohort processing as frozen in docs\ANALYSIS_FREEZE_v1.md.
rem Processes every participant whose memory recording is on disk; finished participants are skipped,
rem so it can be re-run after an interruption. Several hours in total.
rem Needs about 10 GB of free memory for the longest recordings.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -c "import h5py" 2>nul
if errorlevel 1 goto :needh5
"%PY%" -W ignore -m hb.cohort %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo STEP 08 DONE for all participants on disk.
if not "%RC%"=="0" echo STEP 08 finished with failures, exit code %RC%. Send back the newest logs\cohort_*.log
popd
exit /b %RC%
:needh5
echo Twelve recordings above 2 GB use a file format that needs the h5py package.
echo Install it yourself with the line below, then re-run this script:
echo    "%PY%" -m pip install h5py
:fail
popd
exit /b 1
