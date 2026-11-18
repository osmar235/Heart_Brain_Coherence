@echo off
setlocal EnableExtensions
rem Step 13: review C4 and C5, joint pipeline with phase and shift nulls, 16 surrogates each (about 3 to 4 hours).
rem See docs\REVIEW_RESPONSE_PLAN_v1.md. Resumable: finished participants are skipped. Optional: --workers N
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.joint_calibration %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest log in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES, exit code %RC%. Send back the newest log in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
