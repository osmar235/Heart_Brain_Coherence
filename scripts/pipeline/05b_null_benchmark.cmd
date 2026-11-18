@echo off
setlocal EnableExtensions
rem Step 05b: SIMULATION of null model controls on the real trial timing of the pilot participants.
rem Reads only events.tsv files, no recordings. Default 20 simulated sessions per participant and
rem scenario takes roughly 15 minutes. Use --n-sims 4 for a quick look.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.null_benchmark --n-sims 20 %*
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
:fail
popd
exit /b 1
