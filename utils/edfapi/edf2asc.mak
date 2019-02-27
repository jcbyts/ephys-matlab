# Microsoft Developer Studio Generated NMAKE File, Based on edf2asc.dsp
!IF "$(CFG)" == ""
CFG=edf2asc - Win32 Debug
!MESSAGE No configuration specified. Defaulting to edf2asc - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "edf2asc - Win32 Release" && "$(CFG)" != "edf2asc - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "edf2asc.mak" CFG="edf2asc - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "edf2asc - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "edf2asc - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "edf2asc - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\edf2asc.exe" "$(OUTDIR)\edf2asc.bsc"


CLEAN :
	-@erase "$(INTDIR)\AVIEDFTime.obj"
	-@erase "$(INTDIR)\AVIEDFTime.sbr"
	-@erase "$(INTDIR)\main.obj"
	-@erase "$(INTDIR)\main.sbr"
	-@erase "$(INTDIR)\printevent.obj"
	-@erase "$(INTDIR)\printevent.sbr"
	-@erase "$(INTDIR)\printrecordinginfo.obj"
	-@erase "$(INTDIR)\printrecordinginfo.sbr"
	-@erase "$(INTDIR)\printsample.obj"
	-@erase "$(INTDIR)\printsample.sbr"
	-@erase "$(INTDIR)\scenecameraedf.obj"
	-@erase "$(INTDIR)\scenecameraedf.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(OUTDIR)\edf2asc.bsc"
	-@erase "$(OUTDIR)\edf2asc.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\edf2asc.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\edf2asc.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\AVIEDFTime.sbr" \
	"$(INTDIR)\main.sbr" \
	"$(INTDIR)\printevent.sbr" \
	"$(INTDIR)\printrecordinginfo.sbr" \
	"$(INTDIR)\printsample.sbr" \
	"$(INTDIR)\scenecameraedf.sbr"

"$(OUTDIR)\edf2asc.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=..\lib\win32\edfapi.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib setargv.obj /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\edf2asc.pdb" /machine:I386 /out:"$(OUTDIR)\edf2asc.exe" 
LINK32_OBJS= \
	"$(INTDIR)\AVIEDFTime.obj" \
	"$(INTDIR)\main.obj" \
	"$(INTDIR)\printevent.obj" \
	"$(INTDIR)\printrecordinginfo.obj" \
	"$(INTDIR)\printsample.obj" \
	"$(INTDIR)\scenecameraedf.obj"

"$(OUTDIR)\edf2asc.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "edf2asc - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\edf2asc.exe" "$(OUTDIR)\edf2asc.bsc"


CLEAN :
	-@erase "$(INTDIR)\AVIEDFTime.obj"
	-@erase "$(INTDIR)\AVIEDFTime.sbr"
	-@erase "$(INTDIR)\main.obj"
	-@erase "$(INTDIR)\main.sbr"
	-@erase "$(INTDIR)\printevent.obj"
	-@erase "$(INTDIR)\printevent.sbr"
	-@erase "$(INTDIR)\printrecordinginfo.obj"
	-@erase "$(INTDIR)\printrecordinginfo.sbr"
	-@erase "$(INTDIR)\printsample.obj"
	-@erase "$(INTDIR)\printsample.sbr"
	-@erase "$(INTDIR)\scenecameraedf.obj"
	-@erase "$(INTDIR)\scenecameraedf.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(OUTDIR)\edf2asc.bsc"
	-@erase "$(OUTDIR)\edf2asc.exe"
	-@erase "$(OUTDIR)\edf2asc.ilk"
	-@erase "$(OUTDIR)\edf2asc.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "_CONSOLE" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "__NO_JNI__" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\edf2asc.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\edf2asc.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\AVIEDFTime.sbr" \
	"$(INTDIR)\main.sbr" \
	"$(INTDIR)\printevent.sbr" \
	"$(INTDIR)\printrecordinginfo.sbr" \
	"$(INTDIR)\printsample.sbr" \
	"$(INTDIR)\scenecameraedf.sbr"

"$(OUTDIR)\edf2asc.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=..\lib\win32\edfapi.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib setargv.obj /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\edf2asc.pdb" /debug /machine:I386 /out:"$(OUTDIR)\edf2asc.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\AVIEDFTime.obj" \
	"$(INTDIR)\main.obj" \
	"$(INTDIR)\printevent.obj" \
	"$(INTDIR)\printrecordinginfo.obj" \
	"$(INTDIR)\printsample.obj" \
	"$(INTDIR)\scenecameraedf.obj"

"$(OUTDIR)\edf2asc.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

SOURCE="$(InputPath)"
DS_POSTBUILD_DEP=$(INTDIR)\postbld.dep

ALL : $(DS_POSTBUILD_DEP)

# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

$(DS_POSTBUILD_DEP) : "$(OUTDIR)\edf2asc.exe" "$(OUTDIR)\edf2asc.bsc"
   copy Debug\edf2asc.exe ..\..\..\bin
	echo Helper for Post-build step > "$(DS_POSTBUILD_DEP)"

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("edf2asc.dep")
!INCLUDE "edf2asc.dep"
!ELSE 
!MESSAGE Warning: cannot find "edf2asc.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "edf2asc - Win32 Release" || "$(CFG)" == "edf2asc - Win32 Debug"
SOURCE=.\AVIEDFTime.cpp

"$(INTDIR)\AVIEDFTime.obj"	"$(INTDIR)\AVIEDFTime.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\main.c

"$(INTDIR)\main.obj"	"$(INTDIR)\main.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\printevent.c

"$(INTDIR)\printevent.obj"	"$(INTDIR)\printevent.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\printrecordinginfo.c

"$(INTDIR)\printrecordinginfo.obj"	"$(INTDIR)\printrecordinginfo.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\printsample.c

"$(INTDIR)\printsample.obj"	"$(INTDIR)\printsample.sbr" : $(SOURCE) "$(INTDIR)"


SOURCE=.\scenecameraedf.c

"$(INTDIR)\scenecameraedf.obj"	"$(INTDIR)\scenecameraedf.sbr" : $(SOURCE) "$(INTDIR)"



!ENDIF 

