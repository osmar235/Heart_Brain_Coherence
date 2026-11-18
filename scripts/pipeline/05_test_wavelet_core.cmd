@echo off
setlocal EnableExtensions
rem Synthetic benchmarks of the wavelet coherence layer. No real data, no internet.
rem Takes well under a minute. Nothing is installed by this script.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
if not exist "logs" mkdir "logs"
"%PY%" -c "import numpy, scipy, matplotlib, tqdm, pycwt; print('numpy', numpy.__version__, '| scipy', scipy.__version__, '| pycwt', pycwt.__version__)"
if errorlevel 1 goto :missing
"%PY%" -W ignore -m unittest discover -s code\tests -p "test_wavelet_core.py" -v > "logs\test_wavelet_core_latest.txt" 2>&1
set "RC=%ERRORLEVEL%"
type "logs\test_wavelet_core_latest.txt"
echo.
if "%RC%"=="0" echo ALL WAVELET BENCHMARKS PASSED. Output saved in logs\test_wavelet_core_latest.txt
if not "%RC%"=="0" echo SOME BENCHMARKS FAILED. Send back logs\test_wavelet_core_latest.txt
popd
exit /b %RC%
:missing
echo.
echo A required package is missing from this Python environment.
echo To install, run this yourself and then re-run this script:
echo    "%PY%" -m pip install -r requirements.txt
:fail
popd
exit /b 1
