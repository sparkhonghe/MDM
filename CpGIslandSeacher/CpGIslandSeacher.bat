@echo off
SET input_path=D:\c
dir %input_path%\*.fa /B >%input_path%\name.txt
For /F %%a in (%input_path%\name.txt) Do cpgTest %%a
pause
@echo off echo hello! Job is done!

