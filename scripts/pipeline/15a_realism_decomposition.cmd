@echo off
setlocal EnableExtensions
rem Step 15a: realism of the synthetic model (fast) then decomposition of the step 14 bias (about 2 hours).
rem See docs\STEP15_PLAN_v1.md. Resumable: finished participants are skipped. Optional: --workers N
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.trial_shuffle realism %*
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto :done
"%PY%" -W ignore -m hb.trial_shuffle decompose %*
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
