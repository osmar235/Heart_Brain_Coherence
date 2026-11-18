@echo off
setlocal EnableExtensions
rem Step 11: checks 1 and 3, coherence and frozen statistics. Needs step 10.
rem See docs\ROBUSTNESS_PLAN_v1.md. Resumable: finished participants are skipped.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.robust_coherence %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest log in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES, exit code %RC%. Send back the newest log in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
