if exist "cppgoslin\parser\KnownGrammars.h" (del -Force "cppgoslin\parser\KnownGrammars.h")
if exist "src\domain\LipidClasses.cpp" (del -Force "src\domain\LipidClasses.cpp")
if exist "cppgoslin\domain\ClassesEnum.h" (del -Force "cppgoslin\domain\ClassesEnum.h")
del -Force "src\domain\*.o"
del -Force "src\parser\*.o"
del -Force "src\tests\*.o"
if exist "libcppGoslin.dll" (del -Force "libcppGoslin.dll")
if exist "libcppGoslin.a" (del -Force "libcppGoslin.a")
if exist "libcppGoslin.so" (del -Force "libcppGoslin.so")
del -Force "*Test"
if exist "writeGrammarsHeader.exe" (del -Force "writeGrammarsHeader.exe")
if exist "writeLipidEnums.exe" (del -Force "writeLipidEnums.exe")