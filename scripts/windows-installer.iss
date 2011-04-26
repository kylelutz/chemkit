[Files]
Source: ..\lib\chemkit.dll; DestDir: {app}\bin
Source: ..\lib\chemkit-graphics.dll; DestDir: {app}\bin
Source: ..\bin\chemkit-builder.exe; DestDir: {app}\bin
Source: ..\share\chemkit\plugins\*.dll; DestDir: {app}\bin\plugins
Source: ..\share\chemkit\plugins\data\*; DestDir: {app}\bin\plugins\data; Flags: recursesubdirs
Source: C:\Qt\4.6.3\bin\QtCore4.dll; DestDir: {app}\bin
Source: C:\Qt\4.6.3\bin\QtGui4.dll; DestDir: {app}\bin
Source: C:\Qt\4.6.3\bin\QtOpenGL4.dll; DestDir: {app}\bin
Source: ..\lib\blas_win32_MT.dll; DestDir: {app}\bin
Source: ..\lib\lapack_win32_MT.dll; DestDir: {app}\bin
Source: ..\lib\chemkit.lib; DestDir: {app}\lib
Source: ..\lib\chemkit-graphics.lib; DestDir: {app}\lib
Source: ..\src\apps\builder\icons\molecule.ico; DestDir: {app}\icons
Source: ..\src\chemkit\*.h; DestDir: {app}\include\chemkit
Source: ..\src\graphics\*.h; DestDir: {app}\include\chemkit
[Registry]
Root: HKLM; Subkey: SOFTWARE\chemkit; ValueType: string; ValueName: PluginPath; ValueData: {app}\bin\plugins; Flags: uninsdeletekey
Root: HKLM; Subkey: SOFTWARE\chemkit; ValueType: string; ValueName: Version; ValueData: 0.1; Flags: uninsdeletekey
[Setup]
AppName=chemkit
AppVerName=chemkit 0.1
LicenseFile=..\LICENSE
DefaultDirName={pf}\chemkit
DisableProgramGroupPage=true
DefaultGroupName=chemkit
OutputBaseFilename=chemkit-0.1-setup
VersionInfoVersion=0.1
OutputDir=.
[Icons]
Name: {group}\chemkit-builder; Filename: {app}\bin\chemkit-builder.exe; WorkingDir: {app}\bin; IconFilename: {app}\icons\molecule.ico
