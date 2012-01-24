#!/tool/pandora/bin/perl5.8.0

=for

=cut

package neurosim;

# Packages:
use warnings;
use strict;
use DBI;
use POSIX qw(ceil floor);
use POSIX qw(strftime);
use Getopt::Long qw(:config pass_through);
use FindBin;
use Data::Dumper;
use Switch;
use XML::Simple;
use XML::Parser;
use Fcntl;

# Constants:
use constant WIN_COMPILER_VC                                  => "vc";
use constant WIN_COMPILER_MINGW                               => "mingw";
use constant TARGET_ARC_X86_64                                => "x86_64";
use constant OS_WINDOWS                                       => "MSWin32";
use constant OS_LINUX                                         => "linux";

# Globals:
our $CONFIG;
our $TARGET_ARCHITECTURE;
our $OCL_DIR;
our $WIN_COMPILER;
our $WIN_VC_DIR;
our $WIN_VC_ENV_SCRIPT;
our $WIN_MINGW_DIR;
our $WIN_MSYS_DIR;
our $MKDIR = "";
our $COMP = "";
our $LINK = "";
our $RMDIR = "";
our $RM = "";
our $CD = "cd";
our $AR = "";
our $INSTALL = "";
our $I = "";
our $D = "";
our $CURRENT_DIR  = "";
our $CURRENT_REPORT_LOG  = "";
# Find application root directory
our $SCRIPT_ROOT_DIR = $FindBin::Bin; chomp( $SCRIPT_ROOT_DIR );
# Configuration definitions file
our $CONFIG_FILE = "$SCRIPT_ROOT_DIR/config/config.xml";
# Detect OS
our $OS = "$^O\n"; chomp( $OS );

if( not($OS eq OS_WINDOWS) )
{
  die "ERROR: Unsupported OS: $OS\n";
}

# Parse configuration definitions file
if( checkSyntax( $CONFIG_FILE ) ) 
{
  my $xml = new XML::Simple;
  $CONFIG = $xml->XMLin( $CONFIG_FILE );
  #print Dumper($CONFIG);
  
  foreach my $p (keys %{$CONFIG->{config}->{parameters}->{$OS}->{param}})
  {
    my $c = $CONFIG->{config}->{parameters}->{$OS}->{param}->{$p}->{content};
    eval("our \$$p"." = $c;");
  }
}
else 
{
  die "\nError while opening $CONFIG_FILE\n";
}

# Target architecture
if(not($TARGET_ARCHITECTURE eq TARGET_ARC_X86_64))
{
  die "Unsupported target architecture: $TARGET_ARCHITECTURE\n";
  # TODO: when opening 32-bit option need to consider compiler parametes, such as "/MACHINE:X64";
}

# Library build parameters
our $UTIL_NAME = "SDKUtil";
our $UTIL_SRC_DIR = "cpp_cl/util/$UTIL_NAME";
our $UTIL_BIN_DIR = "$UTIL_SRC_DIR/build/release/$TARGET_ARCHITECTURE";
our @UTIL_NAMES = ("SDKApplication", "SDKBitMap", "SDKCommandArgs", "SDKCommon", "SDKFile", 
  "SDKThread");
our @UTIL_INCL_DIR = ("include", "include/$UTIL_NAME", "include/GL");
our $UTIL_INCL_DIR = "";
our $UTIL_COMPILER_OPS = "";
our $UTIL_INSTALL_DIR = "lib/$TARGET_ARCHITECTURE";

# Application build parameters
our $OCL_COMPILER_OPTIONS_FILE = "oclCompilerOptions.txt";
our $APP_NAME = "Neuro";
our $APP_SRC_DIR = "cpp_cl/app/$APP_NAME";
our $APP_BIN_DIR = "$APP_SRC_DIR/bin/release/$TARGET_ARCHITECTURE";
our @APP_INCL_DIR = ("include", "include/$UTIL_NAME", "include/GL");
our $APP_INCL_DIR = "";
our @APP_SRC_NAMES = ("iz_util.c", "integ_util.c", "IntegrationTest.cpp");
our @APP_KERNEL_NAMES = ("Kernel_ExpandEvents.cl", "Kernel_GroupEvents.cl", 
  "Kernel_MakeEventPointers.cl", "Kernel_ScanHistogram.cl", "Kernel_SortSynapticEvents.cl", 
  "Kernel_UpdateNeurons.cl");
our @APP_MISC_FILES = ("Definitions.h", $OCL_COMPILER_OPTIONS_FILE);
our $APP_COMPILER_OPS = "";
our $APP_LINK_OPS = "";
our $APP_LIB_CONFIG = "";
our $APP_INSTALL_DIR = "$SCRIPT_ROOT_DIR/bin/$TARGET_ARCHITECTURE";
our @APP_CURRENT_CONFIG_KEYS = ();
our (%APP_CURRENT_CONFIG); %APP_CURRENT_CONFIG = ();
our (%APP_COMMON_CONFIG); %APP_COMMON_CONFIG = ();
our $APP_COMMON_CONFIG_OPTS = "";
our $APP_REPORT_LOG = "log_report.txt";
our %APP_CURRENT_INFO = ("Test ID"  =>  "",  "SubTest ID"  =>  "");
our %APP_REPORT = ("Compile V"  =>  "",  "Verify"  =>  "", "Compile M"  =>  "", "Run"  =>  "");

# OS-specific builds
if( $OS eq OS_WINDOWS )
{
  # Determine current directory
  $CURRENT_DIR = `$CD`; chomp( $CURRENT_DIR );
  
  # Define OCL root directory
  if( not(-d $OCL_DIR) ) 
  {
    die "Required software OpenCL SDK is not found in $OCL_DIR. ".
      "Install the software and/or modify OCL_DIR constant\n";
  }
  
  # Verify msys path
  if( not(-d $WIN_MSYS_DIR) )
  {
    die "Required software msys is not found in $WIN_MSYS_DIR. ".
      "Install the software and/or modify WIN_MSYS_DIR parameter in config.xml\n";
  }
  else
  {
=for
  TODO:
  Note: To use Msys under MinGW64-x64: 
  a. Open the fstab file (available at msys/1.0/etc/) 
  b. Modify according to your MinGW64-w64 path. For example, modify
  C:\MinGW\ /mingw to C:\MinGW64\ /ming
=cut
  }

  # Set compilation environment based on compiler selection
  my $COMPILER_DIR = "";
  
  if( $WIN_COMPILER eq WIN_COMPILER_MINGW ) 
  {
    if( not(-d $WIN_MINGW_DIR) )
    {
      die "Required software MinGW64 is not found in $WIN_MINGW_DIR. ".
        "Install the software and/or modify WIN_MINGW_DIR parameter in config.xml\n";
    }
    $COMPILER_DIR = $WIN_MINGW_DIR;
    
    $APP_LIB_CONFIG = "-L\"/usr/X11R6/lib\" -lSDKUtil -lOpenCL ".
      "-L\"$SCRIPT_ROOT_DIR/lib/$TARGET_ARCHITECTURE\" ".
      "-L\"$OCL_DIR/lib/$TARGET_ARCHITECTURE\"";
      
    $I = "-I";
    $D = "-D";
    $COMP = "g++";
    $LINK = "g++";
    $UTIL_COMPILER_OPS = "-m64 -Wpointer-arith -Wfloat-equal -g3 -ffor-scope";
    $APP_COMPILER_OPS = "-m64 -Wpointer-arith -Wfloat-equal -g3 -ffor-scope";
    #$APP_COMPILER_OPS = "-m32 -msse2 -Wpointer-arith   -Wfloat-equal -g3 -ffor-scope";
    $APP_LINK_OPS = "";
    $AR = "ar -rsc";

    # Add required entries to the PATH
    print "Adding $WIN_MSYS_DIR; to the PATH\n";
    $ENV{"PATH"} = "$WIN_MSYS_DIR;".$ENV{"PATH"};
    print "Adding $COMPILER_DIR; to the PATH\n";
    $ENV{"PATH"} = "$COMPILER_DIR;".$ENV{"PATH"};
  }
  
  if( $WIN_COMPILER eq WIN_COMPILER_VC )
  {
    if( not(-d $WIN_VC_DIR) )
    {
      die "Required software Visual Studio is not found in $WIN_VC_DIR. ".
        "Install the software and/or modify WIN_VC_DIR parameter in config.xml\n";
    }
    $COMPILER_DIR = $WIN_VC_DIR;
    
    $I = "/I";
    $D = "/D";
    $AR = "\"$COMPILER_DIR/lib.exe\"";
    $COMP = "\"$COMPILER_DIR/cl.exe\"";
    $LINK = "\"$COMPILER_DIR/link.exe\"";
    $UTIL_COMPILER_OPS = "".
      "/c ".            # Compile without linking
      "/Zi ".           # Generate complete debugging information
      "/W2 ".           # Set warning level: 0 - disable; 1-4 (level of severity, high to low); all
      "/WX- ".          # Treat all warnings as errors
      "/O2 ".           # Create fast code
      "/Oi ".           # Replace function calls with intrinsic or special forms for faster code.
      "/GL ".           # Enable whole program optimization
      "/Gm- ".          # Enable minimal rebuild based on .idb file
      "/EHsc ".         # Set exception handling to a - asynchronouse, s - C++ only,
      "/GS ".           # Detect buffer overruns
      "/Gy ".           # Package individual functions in the form of packaged functions (COMDATs).
      "/fp:precise ".   # Set floating-point behavior to precise
      "/Zc:wchar_t ".   # Make wchar_t a native type
      "/Zc:forScope ".  # Create a separate scope for each "for" loop
      "/Gd ".           # Enforce __cdecl function calling convention where possible
      "/MD ".           # Create a multithreaded DLL using MSVCRT.lib
      "/TP ".           # Treat all files named on the command line as C++ source files
      "/errorReport:none ". # Set error reports sent to MicroSoft
      "/D WIN32 /D NDEBUG /D _LIB";
    $APP_COMPILER_OPS = "".
      "/c ".            # Compile without linking
      "/Zi ".           # Generate complete debugging information
      "/W3 ".           # Set warning level: 0 - disable; 1-4 (level of severity, high to low); all
      "/WX- ".          # Treat all warnings as errors
      "/O2 ".           # Create fast code
      "/Oi ".           # Replace function calls with intrinsic or special forms for faster code.
      "/GL ".           # Enable whole program optimization
      "/Gm- ".          # Enable minimal rebuild based on .idb file
      "/EHsc ".         # Set exception handling to a - asynchronouse, s - C++ only,
      "/GS ".           # Detect buffer overruns
      "/Gy ".           # Package individual functions in the form of packaged functions (COMDATs).
      "/fp:precise ".   # Set floating-point behavior to precise
      "/Zc:wchar_t ".   # Make wchar_t a native type
      "/Zc:forScope ".  # Create a separate scope for each "for" loop
      "/Gd ".           # Enforce __cdecl function calling convention where possible
      "/MD ".           # Create a multithreaded DLL using MSVCRT.lib
      "/TP ".           # Treat all files named on the command line as C++ source files
      "/errorReport:none ". # Set error reports sent to MicroSoft
      "/D WIN32 /D NDEBUG /D _CONSOLE /D ATI_OS_WIN /D _CRT_SECURE_NO_DEPRECATE ".
      "/D _CRT_NONSTDC_NO_DEPRECATE";
    $APP_LINK_OPS = "".
      "/errorReport:none ". # Set error reports sent to MicroSoft
      "/INCREMENTAL:NO ".
      "/MANIFEST ".
      "/MANIFESTUAC:\"level='asInvoker' uiAccess='false'\" ".
      "/DEBUG ".
      "/SUBSYSTEM:CONSOLE /OPT:REF /OPT:ICF /LTCG /TLBID:1 /DYNAMICBASE /NXCOMPAT ".
      "/MACHINE:X64";
      
    $APP_LIB_CONFIG = "".
      "/IMPLIB:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.lib\" ".
      "/LIBPATH:\"$SCRIPT_ROOT_DIR/lib/$TARGET_ARCHITECTURE\" ".
      "/LIBPATH:\"$OCL_DIR/lib/$TARGET_ARCHITECTURE\" ".
      "OpenCL.lib SDKUtil.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib ".
      "advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib ".
      "kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ".
      "ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib";
    
    # Set script compilation environment based on $WIN_VC_ENV_SCRIPT
    my $e = `"$WIN_VC_DIR/$WIN_VC_ENV_SCRIPT" & set`; chomp( $e );
    my @env = split("\n", $e);
    foreach $e (@env)
    {
      my @ea = split("=", $e, 2);
      if(scalar(@ea) != 2){next;}
      handleString($ea[0]);
      handleString($ea[1]);
      $ENV{$ea[0]} = $ea[1];
    }
    
    # Add required entries to the PATH
    print "Adding $WIN_MSYS_DIR; to the PATH\n";
    $ENV{"PATH"} = "$WIN_MSYS_DIR;".$ENV{"PATH"};
  }

  # Set additional environment variables
  print "Setting environment variables\n";
  foreach my $p (keys %{$CONFIG->{config}->{environment}->{param}})
  {
    my $c = $CONFIG->{config}->{environment}->{param}->{$p}->{content};
    print "set \"$p=$c\"\n";
    eval("\$ENV{$p} = $c;");
  }
  
  #foreach my $key (sort keys(%ENV)) {print "$key = $ENV{$key}\n";}
  
  # Set compilation tool executable names
  $MKDIR = "mkdir -p";
  $RMDIR = "rm -dfr";
  $RM = "rm";
  $INSTALL = "install -D";
  
  # Create include directive to compiler
  $UTIL_INCL_DIR = "$I \"$OCL_DIR/include\"";
  for my $n (@UTIL_INCL_DIR) 
  {
    $UTIL_INCL_DIR .= " $I \"$SCRIPT_ROOT_DIR/$n\"";
  }
  
  # Create include directive to compiler
  $APP_INCL_DIR = "$I \"$OCL_DIR/include\"";
  for my $n (@APP_INCL_DIR) 
  {
    $APP_INCL_DIR .= " $I \"$SCRIPT_ROOT_DIR/$n\"";
  }
}
else
{
  die "ERROR: Unsupported OS: $OS\n";
}

$| = 1;

main();



####################################################################################################
# main
# Main execution method
####################################################################################################
sub main 
{
  my $test_file_path = "$SCRIPT_ROOT_DIR/config/test.xml";
  my $tests;
  
  # verify existence of test.xml file
  if( not(-e $test_file_path) ) 
  {
    die "\nUnable to locate $test_file_path\n";
  }
  
  # parse refernces to test files
  if( checkSyntax( $test_file_path ) ) 
  {
    my $xml = new XML::Simple;
    $tests = $xml->XMLin( $test_file_path );
    #print Dumper($tests);
  }
  else 
  {
    die "\nError while opening $test_file_path\n";
  }
  
  &compileUtil;

  # sort tests numerically and iterate:
  foreach my $test_id (sort {$a<=>$b} keys %{$tests->{tests}->{test}})
  {
    my $test = $tests->{tests}->{test}->{$test_id};
    
    #skip if marked for skipping 
    if( ($test->{id} eq "skip") || ($test->{path} eq "skip") )
    {
      print("Skipping test reference: id = $test->{id}, path = $test->{path}\n");
      next;
    }
    
    my $config_file_path = "$SCRIPT_ROOT_DIR/".$test->{path};
    
    # verify existence of test file
    if( not(-e $config_file_path) ) 
    {
      print("\nUnable to locate $config_file_path. Skipping.\n");
      next;
    }
  
    my $time_stamp = POSIX::strftime("%Y_%m_%d_%H_%M_%S", localtime);

    my $report_file_name = "report_".$test->{id}."_".$time_stamp.".csv";
    $CURRENT_REPORT_LOG = "$SCRIPT_ROOT_DIR/log/$report_file_name";
    my $config;
    
    # parse current test:
    if( checkSyntax( $config_file_path ) ) 
    {
      my $xml = new XML::Simple;
      $config = $xml->XMLin( $config_file_path );
      #print Dumper($config);
    }
    else 
    {
      die "\nXML error\n";
    }

    # parse test configation common for each iteration
    my $detectedReportLog = 0;
    foreach my $k (sort {lc($a) cmp lc($b)} keys %{$config->{common_config}->{param}})
    {
      my $v = $config->{common_config}->{param}->{$k}->{content};
      $v = eval($v);
      $APP_COMMON_CONFIG_OPTS .= " ".$D." ".$k."=".$v;
      $APP_COMMON_CONFIG{$k} = $v;
    }
    
    # Insert test-specific report keys
    my @temp = keys %{$config->{report_config}->{report_tags}->{param}};
    @APP_REPORT{@temp} = ("") x @temp;
      
    # sort executables numerically and iterate:
    foreach my $exe_name (sort {$a<=>$b} keys %{$config->{specific_config}->{exe}})
    {
      my $exe = $config->{specific_config}->{exe}->{$exe_name};
      
      #skip if marked for skipping 
      if( ($exe->{tag} eq "skip") )
      {
        print("Skipping execution configuration: id = $exe->{id}\n");
        next;
      }
      #print Dumper($exe);
      
      # Parse current configuration
      my $hash = $exe->{config}->{param};
      #print Dumper($hash);
      @APP_CURRENT_CONFIG_KEYS = sort(keys %{$hash});
      @APP_CURRENT_CONFIG{@APP_CURRENT_CONFIG_KEYS} = (0) x @APP_CURRENT_CONFIG_KEYS;
      
      $APP_CURRENT_INFO{"Test ID"} = $test->{id};
      $APP_CURRENT_INFO{"SubTest ID"} = $exe->{id};

      # Create a header in the report file
      logMsg((join(',', "Timestamp", 
        (sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO), 
        (sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG), 
        (sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG),
        (sort {lc($a) cmp lc($b)} keys %APP_REPORT))));

      # Execute configuration sequence recursively
      if(scalar(@APP_CURRENT_CONFIG_KEYS) > 0)
      {
        &compileAndRunAppRecursively(0, $hash);
      }
    }
  }
}



####################################################################################################
# compileAndRunAppRecursively
# Compiles and runs app recursively according to subtest configuration
####################################################################################################
sub compileAndRunAppRecursively
{
  my $p = $_[0];
  my $hash = $_[1];
  
  if($p == scalar(@APP_CURRENT_CONFIG_KEYS))
  {
    return &compileAndRunApp;
  }
  else
  {
    foreach my $v (split(',',$hash->{$APP_CURRENT_CONFIG_KEYS[$p]}->{content}))
    {
      $APP_CURRENT_CONFIG{$APP_CURRENT_CONFIG_KEYS[$p]} = eval($v);
      &compileAndRunAppRecursively(($p+1), $hash);
    }
  }
}



####################################################################################################
# compileAndRunApp
# Compiles and runs app according to current combination of parameters from subtest
####################################################################################################
sub compileAndRunApp
{
  my $config = "";
  my $ts = POSIX::strftime("%Y_%m_%d_%H_%M_%S", localtime);
  
  foreach my $k (@APP_CURRENT_CONFIG_KEYS)
  {
    $config .= (" ".$D." ".$k."=".$APP_CURRENT_CONFIG{$k});
  }

  print "Executing test config: ".$config."\n"; 

  # reset report
  @APP_REPORT{(keys %APP_REPORT)} = ("") x (keys %APP_REPORT);
  
  # Make a 1st run in verification mode
  if(!&compileApp($config." ".$D." SIMULATION_MODE=0 ".$APP_COMMON_CONFIG_OPTS))
  {
    # mark as failed to compile
    $APP_REPORT{"Compile V"} = "FAIL";
    # get values sorted by keys
    my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                  @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                  @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                  @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
    # log message
    logMsg(join(',', $ts, @result));
    return 0;
  }
  else
  {
    # mark as passed to compile
    $APP_REPORT{"Compile V"} = "PASS";
    &runApp;
    
    if(&parseResults)
    {
      # mark as passed to verify
      $APP_REPORT{Verify} = "PASS";
    }
    else
    {
      # mark as failed to verify
      $APP_REPORT{"Verify"} = "FAIL";
      # get values sorted by keys
      my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                    @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                    @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                    @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
      # log message
      logMsg(join(',', $ts, @result));
      return 0;
    }
  }

  # Compile and run the same config for measurements
  if(!&compileApp($config." ".$D." SIMULATION_MODE=4 ".$APP_COMMON_CONFIG_OPTS))
  {
    # mark as failed to compile
    $APP_REPORT{"Compile M"} = "FAIL";
    # get values sorted by keys
    my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                  @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                  @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                  @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
    # log message
    logMsg(join(',', $ts, @result));
    return 0;
  }
  else
  {
    # mark as passed to compile
    $APP_REPORT{"Compile M"} = "PASS";
    &runApp;
    
    if(&parseResults)
    {
      # mark as passed to verify
      $APP_REPORT{"Run"} = "PASS";
      # get values sorted by keys
      my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                    @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                    @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                    @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
      # log message
      logMsg(join(',', $ts, @result));
      return 1;
    }
    else
    {
      # mark as failed to run
      $APP_REPORT{"Run"} = "FAIL";
      # get values sorted by keys
      my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                    @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                    @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                    @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
      # log message
      logMsg(join(',', $ts, @result));
      return 0;
    }
  }
}



####################################################################################################
# compileApp
# Compiles app
####################################################################################################
sub compileApp
{
  my $config = shift;
  my $temp = "";
  
  # Clean
  system( "del /Q \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/\"" );
  system( "del /Q \"$APP_INSTALL_DIR/\"" );
  system( "$RMDIR \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR\"" );
  system( "$MKDIR \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR\"" );

  if( $WIN_COMPILER eq WIN_COMPILER_MINGW ) 
  {
    # Build
    for my $n (@APP_SRC_NAMES) 
    {
      $temp = "$COMP $APP_COMPILER_OPS $config $APP_INCL_DIR ".
        "-o \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$n.o\" ".
        "-c \"$SCRIPT_ROOT_DIR/$APP_SRC_DIR/$n\"";
      print "Building $n\n";
      print "$temp\n";
      system( $temp );
    }

    # Link
    $temp = "$LINK $APP_LINK_OPS -o \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME\"";
    for my $n (@APP_SRC_NAMES) 
    {
      $temp .= " \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$n.o\"";
    }
    $temp .= " $APP_LIB_CONFIG";
    print "Linking\n";
    print "$temp\n";
    system( $temp );
  }
  
  if( $WIN_COMPILER eq WIN_COMPILER_VC ) 
  {
    my $tempS = ""; 
    my $tempO = "";
    for my $n (@APP_SRC_NAMES) 
    {
      $tempS .= "\"$SCRIPT_ROOT_DIR/$APP_SRC_DIR/$n\" ";
      my @pureName = split(/\./, $n);
      $tempO .= "\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/".$pureName[0].".obj\" ";
    }
    
    # Build
    $temp = "$COMP $APP_COMPILER_OPS $config $APP_INCL_DIR ".
      "/Fo\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/\" ".
      "/Fd\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/vc100.pdb\" ".
      "$tempS";
    print "Building application\n";
    print "$temp\n";
    system( $temp );
    
    $temp = "$LINK $APP_LINK_OPS $APP_LIB_CONFIG ".
      "/OUT:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.exe\" ".
      "/ManifestFile:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.exe.intermediate.manifest\" ".
      "/PDB:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/vc100.pdb\" ".
      "$tempO";
    print "Linking application\n";
    print "$temp\n";
    system( $temp );
  }
  
  # Install the app
  print "Installing application $APP_NAME.exe\n";
  system( "$MKDIR \"$APP_INSTALL_DIR\"" );
  $temp = 
    "$INSTALL \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.exe\" \"$APP_INSTALL_DIR/$APP_NAME.exe\"";
  print "$temp\n";
  system( $temp );
  if( not(-e "$APP_INSTALL_DIR/$APP_NAME.exe") ) 
  {
    print "The application was not installed in \"$APP_INSTALL_DIR/$APP_NAME.exe\"\n";
    return 0;
  }
  for my $n ((@APP_MISC_FILES, @APP_KERNEL_NAMES)) 
  {
    $temp = "$INSTALL \"$SCRIPT_ROOT_DIR/$APP_SRC_DIR/$n\" \"$APP_INSTALL_DIR/$n\"";
    print "$temp\n";
    system( $temp );
    if( not(-e "$APP_INSTALL_DIR/$n") ) 
    {
      print "The file $n was not installed in \"$APP_INSTALL_DIR/$n\"\n";
      return 0;
    }
  }
  
  replaceSubStr($config, "/D", "-D");
  my $optionsFile = "$APP_INSTALL_DIR/$OCL_COMPILER_OPTIONS_FILE";
  sysopen(MYOUTFILE, $optionsFile, O_WRONLY|O_APPEND|O_CREAT) || 
    die("compileApp: Cannot open file $optionsFile");
  print MYOUTFILE $config."\n";
  close(MYOUTFILE);
  
  return 1;
}



####################################################################################################
# runApp
# Runs app
####################################################################################################
sub runApp
{
  print "Running $APP_NAME.exe\n";
  chdir("$APP_INSTALL_DIR");
  system("$APP_INSTALL_DIR/$APP_NAME.exe");
  chdir("$APP_INSTALL_DIR");
}



####################################################################################################
# parseResults
# Parses results passed from app  in a report log file
####################################################################################################
sub parseResults
{
  my $key; my $value; my $result = 0;
  
  if( not(-e "$APP_INSTALL_DIR/$APP_REPORT_LOG") ) 
  {
    return $result;
  }

  sysopen(MYINPUTFILE, "$APP_INSTALL_DIR/$APP_REPORT_LOG", O_RDONLY) || 
    die("parseResults: cannot open file $APP_INSTALL_DIR/$APP_REPORT_LOG");
  my(@lines) = <MYINPUTFILE>; 
  close(MYINPUTFILE);
  
  foreach my $line (@lines)
  {
    chomp $line;
    ($key, $value) = split(/:/, $line);

    # detect result
    if($key =~ /^Result/)
    {
      if($value =~ /^PASS/)
      {
        $result = 1;
      }
    }
    
    #allow only predefined keys
    next unless (exists($APP_REPORT{$key}));
    $APP_REPORT{$key} = $value;
  }
  
  return $result;
}



####################################################################################################
# compileUtil
# Compiles utility (library)
####################################################################################################
sub compileUtil
{
  my $temp = "";
  
  # Clean bin dir
  system( "$RMDIR \"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR\"" );
  system( "$MKDIR \"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR\"" );

  if( $WIN_COMPILER eq WIN_COMPILER_MINGW ) 
  {
    # Build libraries
    $temp = "";
    for my $n (@UTIL_NAMES) 
    {
      my $compile = 
        "$COMP $UTIL_COMPILER_OPS $UTIL_INCL_DIR -o \"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/$n.o\" ".
        "-c \"$SCRIPT_ROOT_DIR/$UTIL_SRC_DIR/$n.cpp\"";
      print "Building $n\n";
      print "$compile\n";
      system( $compile );
      $temp .= "\"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/$n.o\" ";
    }
    
    # Create a combined library archive
    print "Creating library lib$UTIL_NAME.a\n";
    $temp = "$AR \"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/lib$UTIL_NAME.a\" $temp";
    print "$temp\n";
    system( $temp );
    
    # Install the library
    print "Installing library lib$UTIL_NAME.a\n";
    $temp = "$INSTALL \"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/lib$UTIL_NAME.a\" ".
      "\"$SCRIPT_ROOT_DIR/$UTIL_INSTALL_DIR/lib$UTIL_NAME.a\"";
    print "$temp\n";
    system( "$RM \"$SCRIPT_ROOT_DIR/$UTIL_INSTALL_DIR/lib$UTIL_NAME.a\"" );
    system( $temp );
  }
  
  if( $WIN_COMPILER eq WIN_COMPILER_VC ) 
  {
    # Build libraries
    my $tempS = ""; 
    my $tempO = "";
    for my $n (@UTIL_NAMES) 
    {
      $tempS .= "\"$SCRIPT_ROOT_DIR/$UTIL_SRC_DIR/$n.cpp\" ";
      $tempO .= "\"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/$n.obj\" ";
    }
    my $compile = "$COMP $UTIL_COMPILER_OPS $UTIL_INCL_DIR ".
      "/Fo\"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/\" ".
      "/Fd\"$SCRIPT_ROOT_DIR/$UTIL_BIN_DIR/vc100.pdb\" ".
      "$tempS";
    print "Building Library\n";
    print "$compile\n";
    system( $compile );

    # Create and install library
    $compile = "$AR /OUT:\"$SCRIPT_ROOT_DIR/$UTIL_INSTALL_DIR/$UTIL_NAME.lib\" /LTCG $tempO";
    print "Creating and installing library $UTIL_NAME.lib\n";
    print "$compile\n";
    system( $compile );
  }
}



####################################################################################################
# checkSyntax
# Checks syntax of an XML file
# 0 -> Failure / 1 -> Success
####################################################################################################
sub checkSyntax
{
	my $file_loc = shift;
	my $xml_parser = XML::Parser->new(ErrorContext=>2, Style=>'Tree');
	# to check well-formedness
	eval { $xml_parser->parsefile($file_loc); };

	if ($@) 
  {
		$@ =~ s/at \/.*?$//s;
		die "checkSyntax(): in $file_loc :\n$@\n";
		# failure
		return 0;
	}
	else 
  {
		## success
		return 1;
	}
}



####################################################################################################
# logMsg
# Logs message into a file
####################################################################################################
sub logMsg
{
	my $msg = $_[0];

  sysopen(MYOUTFILE, $CURRENT_REPORT_LOG, O_WRONLY|O_APPEND|O_CREAT) || 
    die("logMsg: cannot open file $CURRENT_REPORT_LOG");
  print MYOUTFILE $msg."\n";
  close(MYOUTFILE);
}



####################################################################################################
# parseFile
# Parses file into a tree structure
####################################################################################################
sub parseFile
{
	my $file_loc = shift;
	my $simple_parser = XML::Simple->new();
	return $simple_parser->XMLin($file_loc);
}



####################################################################################################
# replaceSubStr
# Replaces a substring of a string with a new substring
####################################################################################################
sub replaceSubStr
{
	$_[0] =~ s/\Q$_[1]\E/$_[2]/g;
}



####################################################################################################
# handleString
# Removes unwanted chars from both ends of string.
####################################################################################################
sub handleString
{
  chomp($_[0]); #remove newline
  $_[0] =~ s/\r//g; #remove carriage return
  $_[0] =~ s/^\s+//; #remove leading spaces
  $_[0] =~ s/\s+$//; #remove trailing spaces
}
