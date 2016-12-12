<?xml version="1.0" encoding="UTF-8" ?>
<config
    Name="MinGW64 Compiler (C)"
    ShortName = "mingw64"
    Manufacturer="GNU"
    Version="4.x"
    Language="C"
    Priority="E"
    Location="$MINGWROOT" >
    <Details
        CompilerExecutable="$CC"
        CompilerDefines="$DEFINES"
        CompilerFlags="$CFLAGS"
        OptimizationFlags="$COPTIMFLAGS"
        DebugFlags="$CDEBUGFLAGS"
        IncludeFlags="$INCLUDE"
        LinkerExecutable="$LD"
        LinkerFlags="$LDFLAGS"
        LinkerLibraries="$LINKLIBS"
        LinkerOptimizationFlags="$LDOPTIMFLAGS"
        LinkerDebugFlags="$LDDEBUGFLAGS"
        CommandLineShell="set MINGW_ROOT_PATH=$MINGWROOT"
        CommandLineShellArg=""
        CompilerDefineFormatter="-D%s"
        LinkerLibrarySwitchFormatter="-l%s;-llib%s"
        LinkerPathFormatter="-L%s"
        LibrarySearchPath="$$LIB;$$LIBPATH;$$PATH;$$INCLUDE;$MATLABROOT\extern\lib\$ARCH\mingw64"
    />
    <vars
          CMDLINE1="$CC -c $DEFINES $INCLUDE $CFLAGS $OPTIM $SRC -o $OBJ"
          CMDLINE2="$LD $LDFLAGS $LDTYPE $LINKOPTIM $LINKEXPORTVER $OBJS $CLIBS $LINKLIBS -o $EXE"
          
          CC="$MINGWROOT\bin\gcc"
          COMPILER="$CC"		  
          DEFINES="-m64 $MATLABMEX"
          MATLABMEX="-DMATLAB_MEX_FILE "
          CFLAGS="-fexceptions -fno-omit-frame-pointer"
          INCLUDE="-I&quot;$MATLABROOT/extern/include&quot; -I&quot;$MATLABROOT/simulink/include&quot; -I&quot;$MATLABROOT/extern\lib\win64\mingw64&quot;"
          COPTIMFLAGS="-O -DNDEBUG"
          CDEBUGFLAGS="-g"
          
          LD="$CC"
		  LINKER="$LD"
          LDTYPE="-shared"
          LDFLAGS="-m64 -Wl,--no-undefined"          
          LINKEXPORT="-Wl,&quot;$MATLABROOT/extern/lib/win64/mingw64/mexFunction.def&quot;"
          LINKEXPORTVER="-Wl,&quot;$MATLABROOT/extern/lib/win64/mingw64/exportsmexfileversion.def&quot;"
          LIBLOC="$MATLABROOT\extern\lib\win64\mingw64"
          LINKLIBS="-L&quot;$MATLABROOT\extern\lib\$ARCH\mingw64&quot; -llibmx -llibmex -llibmat -lm -llibmwlapack -llibmwblas"          
          LDOPTIMFLAGS="-s"
          LDDEBUGFLAGS="-g"

          OBJEXT=".obj"
          LDEXT=".mexw64"              
         
          SETENV="set COMPILER=gcc
				set COMPFLAGS=-c $CFLAGS $DEFINES $MATLABMEX  -fopenmp
				set OPTIMFLAGS=$COPTIMFLAGS 
				set DEBUGFLAGS=$CDEBUGFLAGS 
				set LINKER=gcc 
				set LINKFLAGS=$LDFLAGS $LDTYPE $LINKLIBS $LINKEXPORT 
				set LINKDEBUGFLAGS=$LDDEBUGFLAGS
				set NAME_OUTPUT=-o &quot;%OUTDIR%%MEX_NAME%%MEX_EXT%&quot;"
    />
    <client>
        <engine         
          LINKLIBS="$LINKLIBS -llibeng"
          LINKEXPORT=""
          LINKEXPORTVER=""
          LDEXT=".exe" 
          LINKTYPE=""
          MATLABMEX=""
          LDTYPE=""
        />
    </client>
    <locationFinder>
        <MINGWROOT>
            <and>
                <envVarExists name="MW_MINGW64_LOC" />
                <fileExists name="$$\bin\gcc.exe" />
                <dirExists name="$$\..\" />
            </and>
            
        </MINGWROOT>
    </locationFinder>
    <env
        PATH = "$MINGWROOT\bin;$MATLABROOT\extern\include\$ARCH;$MATLABROOT\extern\include;$MATLABROOT\simulink\include;$MATLABROOT\lib\$ARCH"
        INCLUDE = "$MINGWROOT\include;"
        LIB = "$MINGWROOT\lib;"
        MW_TARGET_ARCH = "win64"        
        LIBPATH="$MATLABROOT\extern\lib\win64"
    />
</config>
