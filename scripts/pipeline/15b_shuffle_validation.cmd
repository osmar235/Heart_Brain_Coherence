@echo off
setlocal EnableExtensions
rem Step 15b: trial shuffle null on synthetic data, independent and matched (about 3 hours). Needs 15a.
rem See docs\STEP15_PLAN_v1.md. Resumable: finished participants are skipped. Optional: --workers N
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.trial_shuffle validate %*
set "RC=%ERRORLEVEL%"
:done
echo.
if "%RC%"=="0" echo DONE. Send back the newest logs in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES OR STOPPED, exit code %RC%. Send back the newest logs in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
