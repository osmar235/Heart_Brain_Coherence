@echo off
setlocal EnableExtensions
rem Downloads EEG and ECG recordings of the pilot participants listed in config\download_config.json, about 4.7 GB.
rem Safe to re-run: finished files are skipped, partial files resume.
rem Extra arguments are passed through, for example --dry-run or --verify-only.
pushd "%~dp0..\.." || exit /b 1
if /i "%HB_CONFIRM_DOWNLOAD%"=="YES" goto :confirmed
echo This step downloads large recordings from OpenNeuro (see the comment at the top of this file).
set /p "ANS=Type YES to start the download: "
if /i not "%ANS%"=="YES" (echo Download not started. & popd & exit /b 4)
:confirmed
call "%~dp0_env.cmd" || goto :fail
"%PY%" -m hb.download_ds008302 --mode pilot %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Log files are in the logs folder.
if not "%RC%"=="0" echo NOT COMPLETE, exit code %RC%. Re-run this script to resume. If it fails again, send back the newest file in logs.
popd
exit /b %RC%
:fail
popd
exit /b 1
