#!/tool/pandora/bin/perl5.8.0

=for
  TODO:
  - release vs debug configurations
  - mingw support
  - compiler warnings cleanup
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
use constant COMPILER_WIN_VC                                  => 0x10000;
use constant COMPILER_WIN_MINGW64                             => 0x10001;
use constant TARGET_ARC_X86_64                                => "x86_64";
use constant OS_WINDOWS                                       => "MSWin32";
use constant OS_LINUX                                         => "linux";
use constant NODE_COMPILE                                     => "Node Compile";
use constant NODE_EXECUTE                                     => "Node Execute";
use constant PREPEND_CUE                                      => "PREPEND CUE";
use constant STD_OUT                                          => "stdout.txt";
use constant STD_ERR                                          => "stderr.txt";

# Globals:
our $CONFIG;
our $TARGET_ARCHITECTURE;
our $OCL_DIR;
our $WIN_VC_DIR;
our $WIN_VC_ENV_SCRIPT;
our $WIN_MINGW_DIR;
our $WIN_MSYS_DIR;
our $MKDIR = "";
our $COMP = "";
our $LINK = "";
our $RMDIR = "";
our $RM = "";
our $CP = "";
our $CD = "cd";
our $AR = "";
our $INSTALL = "";
our $I = "";
our $D = "";
our $CURRENT_DIR  = "";
our $CURRENT_REPORT_LOG  = "";
# Find application root directory
our $SCRIPT_ROOT_DIR = $FindBin::Bin; chomp( $SCRIPT_ROOT_DIR );
# Log directory
our $REPORT_LOG_DIR = "$SCRIPT_ROOT_DIR/log";
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
our $UTIL_BIN_DIR = "$UTIL_SRC_DIR/build/$TARGET_ARCHITECTURE";
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
our $APP_BIN_DIR = "$APP_SRC_DIR/bin/$TARGET_ARCHITECTURE";
our @APP_INCL_DIR = ("include", "include/$UTIL_NAME", "include/GL");
our $APP_INCL_DIR = "";
our @APP_SRC_C = 
  (
    "iz_util.c", 
    "integ_util.c", 
    "Common.cpp", 
    "Data_Connectome.cpp", 
    "Data_SpikeEvents.cpp", 
    "Data_SynapticEvents.cpp", 
    "Operator_Expand.cpp",
    "Operator_Scan.cpp", 
    "Operator_Sort.cpp",
    "Neurosim.cpp"
  );
our @APP_SRC_H = 
  (
    "Definitions.h", 
    "iz_util.h", 
    "integ_util.h", 
    "Common.hpp", 
    "Data_Connectome.hpp",
    "Data_SpikeEvents.hpp", 
    "Data_SynapticEvents.hpp", 
    "Operator_Expand.hpp",
    "Operator_Scan.hpp", 
    "Operator_Sort.hpp", 
    "Neurosim.hpp"
  );
our @APP_KERNEL_FILE_NAMES = 
  (
    "Kernel_ExpandEvents.cl", 
    "Kernel_GroupEvents.cl", 
    "Kernel_MakeEventPointers.cl", 
    "Kernel_ScanHistogram.cl", 
    "Kernel_UpdateNeurons.cl",
    "Kernel_Primitives.cl", 
    "Kernel_Primitives.h"
  );
our $APP_DEFINITION_FILE_NAME = "Definitions.h";
our @APP_CONFIG_FILE_NAMES = ("neuron_variables_sample.csv", $OCL_COMPILER_OPTIONS_FILE);
our $APP_COMPILER_OPS = "";
our $APP_LINK_OPS = "";
our $APP_LIB_CONFIG = "";
our $APP_INSTALL_DIR = "$SCRIPT_ROOT_DIR/bin/$TARGET_ARCHITECTURE";
our @APP_CURRENT_CONFIG_KEYS = ();
our (%APP_CURRENT_CONFIG); %APP_CURRENT_CONFIG = ();
our $APP_FLOW = {};
our (%APP_COMMON_CONFIG); %APP_COMMON_CONFIG = ();
our (%APP_REPORT); %APP_REPORT = ();
our $APP_REPORT_LOG = "log_report.txt";
our %APP_CURRENT_INFO = 
  (
    "Test ID"  =>  "",  
    "SubTest ID"  =>  ""
  );
our %APP_SPECIAL_CONFIG = 
  (
    "COMPILER"  =>  ""
  );

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

  # Set additional environment variables
  #
  print "Setting environment variables\n";
  #
  $ENV{"PATH"} = "$WIN_MSYS_DIR;".$ENV{"PATH"};
  #
  foreach my $p (keys %{$CONFIG->{config}->{environment}->{param}})
  {
    my $c = $CONFIG->{config}->{environment}->{param}->{$p}->{content};
    #print "set \"$p=$c\"\n";
    eval("\$ENV{$p} = $c;");
  }
  
  #foreach my $key (sort keys(%ENV)) {print "$key = $ENV{$key}\n";}
  
  # Set compilation tool executable names
  $MKDIR = "mkdir";
  $RMDIR = "rm -dfr";
  $RM = "rm";
  $CP = "cp";
  $INSTALL = "install -D";
  
  if(not(-d $REPORT_LOG_DIR)) 
  {
    system( "$MKDIR \"$REPORT_LOG_DIR\"" );
    
    if(not(-d $REPORT_LOG_DIR)) 
    {
      die("ERROR: Cannot create directory $REPORT_LOG_DIR");
    }
  }
}
else
{
  die "ERROR: Unsupported OS: $OS\n";
}

# $| is an abbreviation for $OUTPUT_AUTOFLUSH. Setting it to nonzero enables autoflush on the 
# currently-selected file handle, which is STDOUT by default. So the effect is to ensure that 
# print statements and the like output immediately.
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
  
  #&compileUtil;

  # sort tests numerically and iterate:
  foreach my $test_id (sort {$a<=>$b} keys %{$tests->{tests}->{test}})
  {
    my $test = $tests->{tests}->{test}->{$test_id};
    
    #skip if marked for skipping 
    if( ($test->{id} eq "skip") || ($test->{path} eq "skip") )
    {
      #print("Skipping test reference: id = $test->{id}, path = $test->{path}\n");
      next;
    }
    
    print("Executing test flow: $test->{path}\n");
    
    my $config_file_path = "$SCRIPT_ROOT_DIR/".$test->{path};
    
    # verify existence of test file
    if( not(-e $config_file_path) ) 
    {
      print("\nUnable to locate $config_file_path. Skipping.\n");
      next;
    }
  
    my $time_stamp = POSIX::strftime("%Y_%m_%d_%H_%M_%S", localtime);

    my $report_file_name = "report_".$test->{id}."_".$time_stamp.".csv";
    $CURRENT_REPORT_LOG = "$REPORT_LOG_DIR/$report_file_name";
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
    foreach my $k (sort {lc($a) cmp lc($b)} keys %{$config->{common_config}->{param}})
    {
      my $v = eval($config->{common_config}->{param}->{$k}->{content});

      if(!defined($v) || ($v eq ""))
      {
        die "\nDetected undefined value for common configuration parameter: ".$k."\n";
      }
      
      $APP_COMMON_CONFIG{$k} = $v;
    }

    # Insert test-specific report keys
    my @temp = keys %{$config->{report_config}->{report_tags}->{param}};
    @APP_REPORT{@temp} = ("") x @temp;
    
    #print Dumper($config->{specific_config});
      
    # sort executables numerically and iterate:
    foreach my $exe_name (sort {$a<=>$b} keys %{$config->{specific_config}->{exe}})
    {
      my $exe = $config->{specific_config}->{exe}->{$exe_name};
      
      #skip if marked for skipping 
      if( ($exe->{tag} eq "skip") )
      {
        #print("Skipping execution configuration: id = $exe->{id}\n");
        next;
      }
      #print Dumper($exe);
      
      # Parse flow configuration
      $APP_FLOW = $exe->{flow}->{node};
      my @flow_keys = keys %{$APP_FLOW};
      @flow_keys = grep {$_ ne 'Start'} @flow_keys;
      # Correct nodes with single parameters
      foreach my $node (@flow_keys)
      {
        #print ref($APP_FLOW->{$node}->{param});
        if(defined($APP_FLOW->{$node}->{param}->{id}))
        {
          my $id = $APP_FLOW->{$node}->{param}->{id};
          my $content = $APP_FLOW->{$node}->{param}->{content};
          $APP_FLOW->{$node}->{param} = {$id => {'content' => $content}};
        }
      }
      
      # Add node keys to report hash
      @temp = map{NODE_COMPILE.' '.$_} @flow_keys;
      @APP_REPORT{@temp} = ("") x @temp;
      @temp = map{NODE_EXECUTE.' '.$_} @flow_keys;
      @APP_REPORT{@temp} = ("") x @temp;

      # Parse current configuration
      my $config_hash = $exe->{config}->{param};
      #print Dumper($config_hash);
      @APP_CURRENT_CONFIG_KEYS = sort(keys %{$config_hash});
      @APP_CURRENT_CONFIG{@APP_CURRENT_CONFIG_KEYS} = (0) x @APP_CURRENT_CONFIG_KEYS;

      $APP_CURRENT_INFO{"Test ID"} = $test->{id};
      $APP_CURRENT_INFO{"SubTest ID"} = $exe->{id};

      # Create a header in the report file
      my @header = 
        ((sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO), 
        (sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG), 
        (sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG),
        (sort {lc($a) cmp lc($b)} keys %APP_REPORT));
        
      logMsg("\"".(join('","', "Timestamp", @header))."\"");
      
      # Execute configuration sequence recursively
      if(scalar(@APP_CURRENT_CONFIG_KEYS) > 0)
      {
        &compileAndRunAppRecursively(0, $config_hash);
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
    return &executeFlow;
  }
  else
  {
    foreach my $v (split(',',$hash->{$APP_CURRENT_CONFIG_KEYS[$p]}->{content}))
    {
      $v = eval($v);
      
      if(!defined($v) || ($v eq ""))
      {
        die "\nDetected undefined value for flow configuration parameter: ".$p."\n";
      }
    
      $APP_CURRENT_CONFIG{$APP_CURRENT_CONFIG_KEYS[$p]} = $v;
      &compileAndRunAppRecursively(($p+1), $hash);
    }
  }
}



####################################################################################################
# executeFlow
# Executes flow
####################################################################################################
sub executeFlow
{
  my $ts = POSIX::strftime("%Y_%m_%d_%H_%M_%S", localtime);
  my $node = $APP_FLOW->{Start}->{result}->{PASS}->{content};
  my $tempD = "TEMP_D";
  my @special_config_keys = keys %APP_SPECIAL_CONFIG;
  
  # reset report
  @APP_REPORT{(keys %APP_REPORT)} = ("") x (keys %APP_REPORT);
  
  # reset special config
  @APP_SPECIAL_CONFIG{@special_config_keys} = ("") x (keys %APP_SPECIAL_CONFIG);
  
  while(1)
  {
    # Parse param config, resolve confilcts (node overwrites flow overwrites common)
    #
    my $config = "";
    my @node_keys = keys %{$APP_FLOW->{$node}->{param}};
    #
    foreach my $k (@node_keys)
    {
      my $v = eval($APP_FLOW->{$node}->{param}->{$k}->{content});
      
      if(!defined($v) || ($v eq ""))
      {
        die "\nDetected undefined value for node configuration parameter: ".$k."\n";
      }

      if($k ~~ @special_config_keys)
      {
        $APP_SPECIAL_CONFIG{$k} = $v;
      }
      else
      {
        $config .= (" ".$tempD." ".$k."=".$v);
      }
    }
    #
    foreach my $k (@APP_CURRENT_CONFIG_KEYS)
    {
      my $v = $APP_CURRENT_CONFIG{$k};

      if(not($k ~~ @node_keys))
      {
        if($k ~~ @special_config_keys)
        {
          $APP_SPECIAL_CONFIG{$k} = $v;
        }
        else
        {
          $config .= (" ".$tempD." ".$k."=".$v);
        }
      }
    }
    #
    foreach my $k (keys %APP_COMMON_CONFIG)
    {
      my $v = $APP_COMMON_CONFIG{$k};

      if(not($k ~~ @node_keys) && not($k ~~ @APP_CURRENT_CONFIG_KEYS))
      {
        if($k ~~ @special_config_keys)
        {
          $APP_SPECIAL_CONFIG{$k} = $v;
        }
        else
        {
          $config .= (" ".$tempD." ".$k."=".$v);
        }
      }
    }

    # Configure compiler
    &configureCompiler();
    
    # After compiler config $D became set
    $config =~ s/$tempD/$D/g;

    # Execute node
    print "Executing node: ".$node."\n";
    #print "Configuration data:\n".$config."\n";
    
    my $result = &compileAndRunApp($config, $node, $ts);
    
    # Get next node based on results
    if($result == 1)
    {
      $node = $APP_FLOW->{$node}->{result}->{PASS}->{content};
    }
    elsif($result == 0)
    {
      print "Failed to execute node: ".$node."\n";
      $node = $APP_FLOW->{$node}->{result}->{FAIL}->{content};
    }
    else
    {
      # get values sorted by keys
      my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                    @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                    @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                    @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
      # log message
      logMsg("\"".(join('","', $ts, @result))."\"");
      
      die "\nUnexpected result: ".$result."\n";
    }
    
    if($node eq 'Exit')
    {
      # get values sorted by keys
      my @result = (@APP_CURRENT_INFO{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_INFO}, 
                    @APP_COMMON_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_COMMON_CONFIG}, 
                    @APP_CURRENT_CONFIG{sort {lc($a) cmp lc($b)} keys %APP_CURRENT_CONFIG}, 
                    @APP_REPORT{sort {lc($a) cmp lc($b)} keys %APP_REPORT});
      # log message
      logMsg("\"".(join('","', $ts, @result))."\"");
      print "Exiting flow\n";
      return $result;
    }
  }
}



####################################################################################################
# configureCompiler
# Configures compiler and associated environment
####################################################################################################
sub configureCompiler
{
  if(!defined($APP_SPECIAL_CONFIG{"COMPILER"}) || ($APP_SPECIAL_CONFIG{"COMPILER"} eq ""))
  {
    die "\nSpecial variable COMPILER is not set\n";
  }

  # Set compilation environment based on compiler selection
  my $COMPILER_DIR = "";
  
  if( $APP_SPECIAL_CONFIG{"COMPILER"} eq COMPILER_WIN_MINGW64 ) 
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
    print "Adding $COMPILER_DIR; to the PATH\n";
    $ENV{"PATH"} = "$COMPILER_DIR;".$ENV{"PATH"};
  }
  
  if( $APP_SPECIAL_CONFIG{"COMPILER"} eq COMPILER_WIN_VC )
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
    
    my $commonOptions = "".
      "/c ".                # Compile without linking
      "/EHa ".              # Exception handling: a,asynchronouse; s,sync (C++ only)
      "/errorReport:none ". # Error reports sent to MicroSoft
      "/fp:strict ".        # Floating-point behavior: fast, precise, strict
      "/fp:except ".        # Floating-point exception model: except, except-
      "/Gr ".               # Function call convention: d,__cdecl; r,__fastcall; z,__stdcall 
      "/GL ".               # Enable whole program optimization (Can't use /Z7, /Zi, /ZI with it)
      "/Gm- ".              # Rebuild based on .idb file: m,enable; m-,disable.
      "/GS ".               # Detect buffer overruns
      "/Gy ".               # Package functions in the form of packaged functions (COMDATs)
      "/MD ".               # Create a multithreaded DLL using MSVCRT.lib
      "/nologo ".           # Suppresses the display of the copyright banner informational messages.
      "/O2 ".               # Optimizations: 1,short code; 2,fast code
      "/TP ".               # Treat all files named on the command line as C++ source files
      #"/Wall ".             # Enable all warnings, including those disabled by default
      "/W3 ".               # Set warning level: 0,disable; 1-4,(level of severity, high to low)
      "/WX ".               # Treat all warnings as errors
      "/wd4710 ".           # Disable warning C4710 (function not inlined)
      "/wd4711 ".           # Disable warning C4711 (selected for automatic inline expansion)
      "/wd4820 ".           # Disable warning C4820 (bytes padding added after data member)
      "/Zc:wchar_t  ".      # wchar_t is a native type: wchar_t,enable; wchar_t-,disable
      "/Zc:forScope  ".     # Create a separate scope for each "for" loop
      #"/Zi ".              # Generate complete debugging information
      " ";
      
    $UTIL_COMPILER_OPS = $commonOptions.
      "/D WIN32 /D NDEBUG /D _LIB /D COMPILER=".COMPILER_WIN_VC;
      
    $APP_COMPILER_OPS = $commonOptions.
      "/D WIN32 /D NDEBUG /D _CONSOLE /D ATI_OS_WIN /D _CRT_SECURE_NO_DEPRECATE ".
      "/D _CRT_NONSTDC_NO_DEPRECATE /D COMPILER=".COMPILER_WIN_VC;
      
    $APP_LINK_OPS = "".
      "/nologo ".       # Suppresses the display of the copyright banner informational messages.
      "/errorReport:none ". # Set error reports sent to MicroSoft
      "/INCREMENTAL:NO ".
      "/MANIFEST ".
      "/MANIFESTUAC:\"level='asInvoker' uiAccess='false'\" ".
      "/DEBUG ".
      "/SUBSYSTEM:CONSOLE /OPT:REF /OPT:ICF /LTCG /TLBID:1 /DYNAMICBASE /NXCOMPAT ".
      "/MACHINE:X64";
      
    $APP_LIB_CONFIG = "".
      "/nologo ".       # Suppresses the display of the copyright banner informational messages.
      "/IMPLIB:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.lib\" ".
      "/LIBPATH:\"$SCRIPT_ROOT_DIR/lib/$TARGET_ARCHITECTURE\" ".
      "/LIBPATH:\"$OCL_DIR/lib/$TARGET_ARCHITECTURE\" ".
      #"OpenCL.lib SDKUtil.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib ".
      "OpenCL.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib ".
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
  }
  
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



####################################################################################################
# compileAndRunApp
# Compiles and runs app according to current combination of parameters from subtest
####################################################################################################
sub compileAndRunApp
{
  my $config = $_[0];
  my $node = $_[1];
  my $ts = $_[2];
  
  # Attempt to compile
  my $compile_error = &compileApp($config);
  #
  if(not($compile_error eq ""))
  {
    # mark as failed to compile
    $APP_REPORT{NODE_COMPILE.' '.$node} = $compile_error;
    return 0;
  }
  else
  {
    # mark as passed to compile
    $APP_REPORT{NODE_COMPILE.' '.$node} = "PASS";
    &runApp;
    
    my $runResult = &parseResults;
    
    if($runResult == 1)
    {
      # mark as passed to execute
      $APP_REPORT{NODE_EXECUTE.' '.$node} = "PASS";
    }
    elsif($runResult == 0)
    {
      # mark as failed to execute with unknown failure cause
      $APP_REPORT{NODE_EXECUTE.' '.$node} = "FAIL: Unknown";
      return 0;
    }
    elsif($runResult == -1)
    {
      # mark as failed to execute with known failure cause
      if(exists($APP_REPORT{"Result"}))
      {
        $APP_REPORT{NODE_EXECUTE.' '.$node} = $APP_REPORT{"Result"};
      }
      else
      {
        $APP_REPORT{NODE_EXECUTE.' '.$node} = "FAIL: Cannot identify cause";
      }
      return 0;
    }
  }

  return 1;
}



####################################################################################################
# compileApp
# Compiles app
####################################################################################################
sub compileApp
{
  my $config = shift;
  my $temp = "";
  my $result = "";
  
  # Clean
  #
  $temp = "del /Q \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/\"";
  $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
  #
  $temp = "del /Q \"$APP_INSTALL_DIR/\"";
  $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
  #
  $temp = "$RMDIR \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR\"";
  $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
  #
  $temp = "$MKDIR \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR\"";
  $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
  
  # Option COMPILER_WIN_MINGW64
  #
  if( $APP_SPECIAL_CONFIG{"COMPILER"} eq COMPILER_WIN_MINGW64 )
  {
    # Build
    #
    for my $n (@APP_SRC_C) 
    {
      $temp = "$COMP $APP_COMPILER_OPS $config $APP_INCL_DIR ".
        "-o \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$n.o\" ".
        "-c \"$SCRIPT_ROOT_DIR/$APP_SRC_DIR/$n\"";
      print "Building $n\n";
      print "$temp\n";
      $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
    }

    # Link
    #
    $temp = "$LINK $APP_LINK_OPS -o \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME\"";
    for my $n (@APP_SRC_C) 
    {
      $temp .= " \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$n.o\"";
    }
    $temp .= " $APP_LIB_CONFIG";
    print "Linking\n";
    print "$temp\n";
    $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
  }
  
  # Option COMPILER_WIN_VC
  #
  if( $APP_SPECIAL_CONFIG{"COMPILER"} eq COMPILER_WIN_VC ) 
  {
    # Copy src files into build dir
    #
    for my $n (@APP_SRC_C, @APP_SRC_H) 
    {
      $temp = "$CP \"$SCRIPT_ROOT_DIR/$APP_SRC_DIR/$n\" \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/\"";
      $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
    }
    
    # Insert configuration-specific defenitions at predetermined cue
    #
    {
      $temp = "";
      
      my @configArray = split("/D", $config);
      
      shift(@configArray); 
      
      for my $n (@configArray) 
      {
        replaceSubStr($n, "=", " ");
        replaceSubStr($n, "\\\"", "\"");
        
        $temp .= "#define $n\n";
      }
      
      my $optionsFile = "$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_DEFINITION_FILE_NAME";
      sysopen(MYINPUTFILE, $optionsFile, O_RDONLY) || 
        die("compileApp: Cannot open file $optionsFile for read");
      my @configFile = <MYINPUTFILE>; 
      close(MYINPUTFILE);
      
      #&validateDefinitions(\@configFile);
      
      my $inserted = 0;
      
      foreach my $i (0..(scalar(@configFile)-1)) 
      {
        if (index($configFile[$i], PREPEND_CUE) != -1)
        {
          $configFile[$i] = "\n".$temp."\n"; $inserted = 1; last;
        }
      }
      
      if($inserted == 0)
      {
        die("compileApp: Cannot insert configuration-specific defenitions in ".
          "$APP_DEFINITION_FILE_NAME");
      }
      
      sysopen(MYOUTFILE, $optionsFile, O_WRONLY|O_TRUNC|O_CREAT) || 
        die("compileApp: Cannot open file $optionsFile for write");
      print MYOUTFILE @configFile;
      close(MYOUTFILE);
    }

    # Assemble src and obj file paths
    #
    my $tempS = ""; my $tempO = "";
    #
    for my $n (@APP_SRC_C) 
    {
      $tempS .= "\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$n\" ";
      my @pureName = split(/\./, $n);
      $tempO .= "\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/".$pureName[0].".obj\" ";
    }
    
    # Build
    #
    $temp = "$COMP $APP_COMPILER_OPS $APP_INCL_DIR ".
      "/Fo\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/\" ".
      "/Fd\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/vc100.pdb\" ".
      "$tempS";
    #
    print "Building application\n";
    #print "$temp\n";
    $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
    
    # Link
    #
    $temp = "$LINK $APP_LINK_OPS $APP_LIB_CONFIG ".
      "/OUT:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.exe\" ".
      "/ManifestFile:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$APP_NAME.exe.intermediate.manifest\" ".
      "/PDB:\"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/vc100.pdb\" ".
      "$tempO";
    #
    print "Linking application\n";
    #print "$temp\n";
    $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
  }
  
  # Install the app
  #
  {
    print "Installing application\n";
    
    if( not(-e "$APP_INSTALL_DIR") ) 
    {
      $temp = "$MKDIR \"$APP_INSTALL_DIR\"";
      $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");
    }
    
    for my $n ("$APP_NAME.exe", $APP_DEFINITION_FILE_NAME)
    {
      $temp = "$INSTALL \"$SCRIPT_ROOT_DIR/$APP_BIN_DIR/$n\" \"$APP_INSTALL_DIR/$n\"";
      #print "$temp\n";
      $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");

      if( not(-e "$APP_INSTALL_DIR/$n") ) 
      {
        print "The file $n was not installed in \"$APP_INSTALL_DIR/$n\"\n";
        return "FAIL: Install 1";
      }
    }
    
    for my $n (@APP_KERNEL_FILE_NAMES)
    {
      $temp = "$INSTALL \"$SCRIPT_ROOT_DIR/$APP_SRC_DIR/$n\" \"$APP_INSTALL_DIR/$n\"";
      #print "$temp\n";
      $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");

      if( not(-e "$APP_INSTALL_DIR/$n") ) 
      {
        print "The file $n was not installed in \"$APP_INSTALL_DIR/$n\"\n";
        return "FAIL: Install 2";
      }
    }
    
    for my $n (@APP_CONFIG_FILE_NAMES) 
    {
      $temp = "$INSTALL \"$SCRIPT_ROOT_DIR/config/$n\" \"$APP_INSTALL_DIR/$n\"";
      #print "$temp\n";
      $result = &systemCall($temp, STD_OUT, STD_ERR); return $result unless($result eq "");

      if( not(-e "$APP_INSTALL_DIR/$n") ) 
      {
        print "The file $n was not installed in \"$APP_INSTALL_DIR/$n\"\n";
        return "FAIL: Install 3";
      }
    }
  }

  return $result;
}



####################################################################################################
# runApp
# Runs app
####################################################################################################
sub runApp
{
  print "Giving execution control to $APP_NAME.exe\n";
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
    #print("Line: $line\n");
    ($key, $value) = split(/:/, $line, 2);

    if(defined($key))
    {
      # detect result
      if($key =~ /^Result/)
      {
        if($value =~ /^PASS/)
        {
          $result = 1;
        }
        if($value =~ /^FAIL/)
        {
          $result = -1;
        }
      }
      
      #allow only predefined keys
      next unless (exists($APP_REPORT{$key}));
      replaceSubStr($value, ",", ";");
      $APP_REPORT{$key} = $value;
    }
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

  if( $APP_SPECIAL_CONFIG{"COMPILER"} eq COMPILER_WIN_MINGW64 ) 
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
  
  if( $APP_SPECIAL_CONFIG{"COMPILER"} eq COMPILER_WIN_VC ) 
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
    #print "$compile\n";
    system( $compile );

    # Create and install library
    $compile = "$AR /OUT:\"$SCRIPT_ROOT_DIR/$UTIL_INSTALL_DIR/$UTIL_NAME.lib\" /LTCG $tempO";
    print "Creating and installing library $UTIL_NAME.lib\n";
    print "$compile\n";
    system( $compile );
  }
}



####################################################################################################
# validateDefinitions
# Validates code in Definitions preprocessor file
####################################################################################################
sub validateDefinitions
{
  my $definitions = $_[0];
  
  foreach my $i (reverse(0..scalar(@$definitions)-1)) 
  {
    my $definition = $definitions->[$i];
    
    # find a #define
    if(index($definition, "#define") != -1)
    {
      # extract definition name for this #define
      chomp( $definition );
      my @linePcs = split(" ", $definition); 
      $definition = $linePcs[1]; chomp( $definition );

      # check if this definition name is used before it is defined by #define
      foreach my $j (reverse(0..($i-1))) 
      {
        my $line = $definitions->[$j];
        
        if(index($line, $definition) != -1)
        {
          chomp( $line );
          next if(index($line, "#if !defined") != -1);
          next if(index($line, "#if !(defined") != -1);
          next if(index($line, "#undef") != -1);
          next if(index($line, "#define") != -1);
          next if(index($line, "#ifndef") != -1);
          next if(index($line, "\/\*") != -1);
          next if(index($line, "\*\/") != -1);
          print $definition." \@ ".($i+1)." -> ".$line." \@ ".($j+1)."\n";
        }
      }
    }
  }
  
  die("Debug");
}



####################################################################################################
# systemCall
# Wrapper around system
####################################################################################################
sub systemCall
{
  my $command = $_[0];
  my $stdout = $_[1];
  my $stderr = $_[2];
  
  my $result = 0;
  
  if(not($stdout eq "") || not($stderr eq ""))
  {
    if(not(-d $REPORT_LOG_DIR)) 
    {
      my $r = 0;
      $r = system( "$MKDIR \"$REPORT_LOG_DIR\"" );
      
      if(not(-d $REPORT_LOG_DIR) || ($r != 0)) 
      {
        die("systemCall: Cannot create directory $REPORT_LOG_DIR. Exit code $r");
      }
    }
  }
  
  if(not($stdout eq ""))
  {
    $stdout = "$REPORT_LOG_DIR/$stdout";
    $command .= " 1>\"$stdout\"";
  }
  
  if(not($stderr eq ""))
  {
    $stderr = "$REPORT_LOG_DIR/$stderr";
    $command .= " 2>\"$stderr\"";
  }
 
  #print $command."\n";
  
  $result = system($command);

	if($result != 0)
  {
    $result = "FAIL: Exist status $result";

    for my $fp ($stdout, $stderr)
    {
      if(not($fp eq ""))
      {
        sysopen(MYINPUTFILE, $fp, O_RDONLY) || die("systemCall: Cannot open file $fp");
        my @file = <MYINPUTFILE>; 
        close(MYINPUTFILE);
        
        $result .= "\n";
        
        foreach my $line (@file)
        {
          $line =~ tr/,/;/;
          $result .= "$line";
        }
      }
    }

		return $result;
	}
  
  return "";
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
