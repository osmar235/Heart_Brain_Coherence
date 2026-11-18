@echo off
setlocal EnableExtensions
rem Steps 16b then 16c in one go (about 8 hours). Stops if 16b has failures.
call "%~dp016b_ica_apply.cmd" || exit /b 1
call "%~dp016c_ica_analyse.cmd"
exit /b %ERRORLEVEL%
