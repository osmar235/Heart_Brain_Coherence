@echo off
setlocal EnableExtensions
rem Step 18, synthetic part: cross participant null validated on synthetic data without coupling, and the
rem positive control (planted coupling at three levels) with the trial shuffle and cross participant nulls.
rem About 9 hours with 6 workers. See docs\STEP18_PLAN_v1.md. Resumable: finished participants are skipped.
rem Optional: --workers N (for example --workers 12 on a machine with 16 cores).
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.symmetric_tests synthetic %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest logs in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES OR STOPPED, exit code %RC%. Send back the newest logs in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
