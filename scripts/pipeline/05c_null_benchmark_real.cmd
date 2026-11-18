@echo off
setlocal EnableExtensions
rem Step 05c: null model benchmark rebuilt on the pilot's own response shapes and noise.
rem Needs step 06 outputs (derivatives\slow). Builds SURROGATE sessions with no ongoing coupling by
rem construction; the real heart and brain coherence is not computed. Roughly 15 to 20 minutes.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -W ignore -m hb.null_benchmark_real %*
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
:fail
popd
exit /b 1
