@echo off
setlocal EnableExtensions
rem Downloads all small BIDS text files of ds008302 v1.0.0, about 11 MB, no recordings.
rem Safe to re-run: finished files are skipped, partial files resume.
rem Extra arguments are passed through, for example --dry-run or --verify-only.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -m hb.download_ds008302 --mode metadata %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Log files are in the logs folder.
if not "%RC%"=="0" echo NOT COMPLETE, exit code %RC%. Re-run this script to resume. If it fails again, send back the newest file in logs.
popd
exit /b %RC%
:fail
popd
exit /b 1
