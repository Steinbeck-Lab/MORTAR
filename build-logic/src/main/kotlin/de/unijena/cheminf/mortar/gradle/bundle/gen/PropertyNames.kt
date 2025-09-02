/*
 * MORTAR - MOlecule fRagmenTAtion fRamework
 * Copyright (C) 2025 Felix Baensch, Jonas Schaub (felix.j.baensch@gmail.com, jonas.schaub@uni-jena.de)
 *
 * Source code is available at <https://github.com/FelixBaensch/MORTAR>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package de.unijena.cheminf.mortar.gradle.bundle.gen

import de.unijena.cheminf.mortar.gradle.bundle.PropertyKey

// Auto-generated at build time - do not edit!
// Generated from: MortarBundle.properties
// Generated at: 2025-09-02T09:14:11.889790500

object PropertyNames {
    /** 
     * Key for the app version property.
     * 
     * Property value: 1.4.0.0
     */
    val APP_VERSION = PropertyKey("app.version")

    /** 
     * Key for the icon path linux property.
     * 
     * Property value: Images/MORTAR-icon.icns
     */
    val ICON_PATH_LINUX = PropertyKey("icon.path.linux")

    /** 
     * Key for the icon path mac property.
     * 
     * Property value: Images/MORTAR-icon.icns
     */
    val ICON_PATH_MAC = PropertyKey("icon.path.mac")

    /** 
     * Key for the icon path windows property.
     * 
     * Property value: Images/Mortar_Logo_Icon1.ico
     */
    val ICON_PATH_WINDOWS = PropertyKey("icon.path.windows")

    /** 
     * Key for the local deploy about url property.
     * 
     * Property value: https://github.com/FelixBaensch/MORTAR
     */
    val LOCAL_DEPLOY_ABOUT_URL = PropertyKey("local.deploy.about.url")

    /** 
     * Key for the local deploy copyright property.
     * 
     * Property value: Copyright (C) 2025  Felix Baensch, Jonas Schaub; licensed under MIT license
     */
    val LOCAL_DEPLOY_COPYRIGHT = PropertyKey("local.deploy.copyright")

    /** 
     * Key for the local deploy desc property.
     * 
     * Property value: MORTAR - MOlecule fRagmenTAtion fRamework
     */
    val LOCAL_DEPLOY_DESC = PropertyKey("local.deploy.desc")

    /** 
     * Key for the local deploy gradle group property.
     * 
     * Property value: deployment
     */
    val LOCAL_DEPLOY_GRADLE_GROUP = PropertyKey("local.deploy.gradle.group")

    /** 
     * Key for the local deploy linux desc property.
     * 
     * Property value: Local Linux packaging via jpackage (produces DEB or RPM)
     */
    val LOCAL_DEPLOY_LINUX_DESC = PropertyKey("local.deploy.linux.desc")

    /** 
     * Key for the local deploy mac desc property.
     * 
     * Property value: Local macOS packaging via jpackage (produces DMG)
     */
    val LOCAL_DEPLOY_MAC_DESC = PropertyKey("local.deploy.mac.desc")

    /** 
     * Key for the local deploy unix java option 1 property.
     * 
     * Property value: -Xms512m
     */
    val LOCAL_DEPLOY_UNIX_JAVA_OPTION_1 = PropertyKey("local.deploy.unix.java.option.1")

    /** 
     * Key for the local deploy unix java option 2 property.
     * 
     * Property value: -Xmx4g
     */
    val LOCAL_DEPLOY_UNIX_JAVA_OPTION_2 = PropertyKey("local.deploy.unix.java.option.2")

    /** 
     * Key for the local deploy win app id yaml path property.
     * 
     * Property value: build-logic\src\main\resources\inno_setup\app-ids.yml
     */
    val LOCAL_DEPLOY_WIN_APP_ID_YAML_PATH = PropertyKey("local.deploy.win.app.id.yaml.path")

    /** 
     * Key for the local deploy win desc property.
     * 
     * Property value: Local Windows deployment (Inno Setup 6 required)
     */
    val LOCAL_DEPLOY_WIN_DESC = PropertyKey("local.deploy.win.desc")

    /** 
     * Key for the local deploy win inno setup compiler path property.
     * 
     * Property value: C:\Program Files (x86)\Inno Setup 6\iscc.exe
     */
    val LOCAL_DEPLOY_WIN_INNO_SETUP_COMPILER_PATH = PropertyKey("local.deploy.win.inno.setup.compiler.path")

    /** 
     * Key for the local deploy win jre dir name property.
     * 
     * Property value: jdk-21.0.1_12_jre
     */
    val LOCAL_DEPLOY_WIN_JRE_DIR_NAME = PropertyKey("local.deploy.win.jre.dir.name")

    /** 
     * Key for the local deploy win jre extracted dir name property.
     * 
     * Property value: jdk-21.0.1+12-jre
     */
    val LOCAL_DEPLOY_WIN_JRE_EXTRACTED_DIR_NAME = PropertyKey("local.deploy.win.jre.extracted.dir.name")

    /** 
     * Key for the local deploy win jre url property.
     * 
     * Property value: https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.1%2B12/OpenJDK21U-jre_x64_windows_hotspot_21.0.1_12.zip
     */
    val LOCAL_DEPLOY_WIN_JRE_URL = PropertyKey("local.deploy.win.jre.url")

    /** 
     * Key for the local deploy win resources dir property.
     * 
     * Property value: build-logic/src/main/resources
     */
    val LOCAL_DEPLOY_WIN_RESOURCES_DIR = PropertyKey("local.deploy.win.resources.dir")

    /** 
     * Key for the msg inno setup not found property.
     * 
     * Property value: Inno Setup not found. Please install Inno Setup 6 first.
     */
    val MSG_INNO_SETUP_NOT_FOUND = PropertyKey("msg.inno.setup.not.found")

    /** 
     * Key for the msg jpackage not found property.
     * 
     * Property value: jpackage not found. Ensure a JDK with jpackage is installed and on PATH.
     */
    val MSG_JPACKAGE_NOT_FOUND = PropertyKey("msg.jpackage.not.found")

    /** 
     * Key for the msg local deploy run on linux property.
     * 
     * Property value: This task must be run on Linux.
     */
    val MSG_LOCAL_DEPLOY_RUN_ON_LINUX = PropertyKey("msg.local.deploy.run.on.linux")

    /** 
     * Key for the msg local deploy run on mac property.
     * 
     * Property value: This task must be run on macOS.
     */
    val MSG_LOCAL_DEPLOY_RUN_ON_MAC = PropertyKey("msg.local.deploy.run.on.mac")

    /** 
     * Key for the msg local deploy run on windows property.
     * 
     * Property value: This task must be run on Windows.
     */
    val MSG_LOCAL_DEPLOY_RUN_ON_WINDOWS = PropertyKey("msg.local.deploy.run.on.windows")

    /** 
     * Key for the os linux property.
     * 
     * Property value: Linux
     */
    val OS_LINUX = PropertyKey("os.linux")

    /** 
     * Key for the os mac property.
     * 
     * Property value: Mac
     */
    val OS_MAC = PropertyKey("os.mac")

    /** 
     * Key for the os win property.
     * 
     * Property value: Windows
     */
    val OS_WIN = PropertyKey("os.win")
}
