@echo off
setlocal EnableExtensions
rem Step 16c: trial shuffle analysis of v2, v2_noica and v2_emgadj (about 6 to 7 hours). Needs 16b.
rem See docs\ANALYSIS_FREEZE_v2_EEG.md. Resumable: finished participants are skipped. Optional: --workers N
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.ica_v2_apply analyse %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest log in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES, exit code %RC%. Send back the newest log in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
