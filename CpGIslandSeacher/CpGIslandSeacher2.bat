@echo off
SET input_path=D:\c
dir %input_path%\*.fna /B >%input_path%\name.txt
For /F %%a in (%input_path%\name.txt) Do cpglh_new %%a
pause
@echo off echo hello! Job is done!

