@echo off
setlocal EnableExtensions
rem Step 18, real part: cross participant null on the real data (v1, ICA v2, ICA v2 with muscle proxy),
rem then the injected coupling test on the real ICA v2 data. Needs 18a (its gate must have passed).
rem About 9 hours with 6 workers. See docs\STEP18_PLAN_v1.md. Resumable. Optional: --workers N
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.symmetric_tests real %*
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto :done
"%PY%" -W ignore -m hb.symmetric_tests inject %*
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
