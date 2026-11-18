@echo off
setlocal EnableExtensions
rem Step 12: checks 2 and 4, cohort scale surrogates in parallel (about 1 to 2 hours). Independent of steps 10 and 11.
rem See docs\ROBUSTNESS_PLAN_v1.md. Resumable: finished participants are skipped.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.null_calibration %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest log in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES, exit code %RC%. Send back the newest log in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
