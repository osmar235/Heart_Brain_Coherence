@echo off
setlocal EnableExtensions
rem Step 14: review C1 C2 C5, independent synthetic null cohorts, 2 scenarios (about 2 hours).
rem See docs\REVIEW_RESPONSE_PLAN_v1.md. Resumable: finished participants are skipped. Optional: --workers N
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.null_validation %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest log in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES, exit code %RC%. Send back the newest log in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
