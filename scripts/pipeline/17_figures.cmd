@echo off
setlocal EnableExtensions
rem Step 17: manuscript figures from the result tables (seconds). No recording is read.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.make_figures %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Figures in results\figures.
popd
exit /b %RC%
:fail
popd
exit /b 1
