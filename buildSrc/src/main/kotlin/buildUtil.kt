import java.io.File

/**
 * Modifies the start scripts for the Mortar application to align with specific configurations
 * for both Windows and Unix systems. The method applies various replacements in the script
 * content to update Java runtime settings, classpath definitions, and other executable
 * configurations.
 *
 * @param aWindowsFile the File object representing the Windows start script to be modified
 * @param anUnixFile the File object representing the Unix start script to be modified
 * @param aMortarOptimizationVariable a specific JVM option to be added or replaced in the scripts
 */
fun modifyMortarStartScripts(aWindowsFile: File, anUnixFile: File, aMortarOptimizationVariable: String) {
    // Windows script edits
    aWindowsFile.writeText(
        aWindowsFile.readText()
            .replace("java.exe", "javaw.exe")
            .replace(
                "if defined JAVA_HOME goto findJavaFromJavaHome",
                "goto setJavaFromAppHome\n\n@rem unused because Java home is set in method called above"
            )
            .replace(":findJavaFromJavaHome", ":setJavaFromAppHome")
            .replace(
                "set JAVA_HOME=%JAVA_HOME:\"=%",
                "set JAVA_HOME=%APP_HOME%\\jdk-21.0.1_12_jre\\"
            )
            .replace(Regex("set CLASSPATH=.*"), "set CLASSPATH=.;%APP_HOME%/lib/*")
            .replace(
                "\"%JAVA_EXE%\" %DEFAULT_JVM_OPTS% %JAVA_OPTS% $aMortarOptimizationVariable  -classpath \"%CLASSPATH%\" de.unijena.cheminf.mortar.main.Main %*",
                "start \"MORTAR\" \"%JAVA_EXE%\" %DEFAULT_JVM_OPTS% %JAVA_OPTS% $aMortarOptimizationVariable  -classpath \"%CLASSPATH%\" de.unijena.cheminf.mortar.main.Main \"-skipJavaVersionCheck\""
            )
            .replace(
                "Please set the JAVA_HOME variable in your environment to match the",
                "Please check your MORTAR installation,"
            )
            .replace(
                "location of your Java installation.",
                "something must be wrong with the Java Runtime Environment shipped with MORTAR."
            )
    )

    // Unix script edits
    anUnixFile.writeText(
        anUnixFile.readText()
            .replace(Regex("CLASSPATH=\$APP_HOME/lib.*"), "CLASSPATH=\$APP_HOME/lib/*")
    )
}
