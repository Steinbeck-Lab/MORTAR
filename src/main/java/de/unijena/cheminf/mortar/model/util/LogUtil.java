/*
 * MORTAR - MOlecule fRagmenTAtion fRamework
 * Copyright (C) 2020  Felix Baensch, Jonas Schaub (felix.baensch@w-hs.de, jonas-schaub@uni-jena.de)
 *
 * Source code is available at <https://github.com/FelixBaensch/MORTAR>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.unijena.cheminf.mortar.model.util;

import de.unijena.cheminf.mortar.gui.util.GuiUtil;
import de.unijena.cheminf.mortar.message.Message;
import javafx.scene.control.Alert;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.LogManager;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.SimpleFormatter;

/**
 * Logging utilities. The Java-own logging API is employed.
 *
 * @author Jonas Schaub
 */
public final class LogUtil {
    //<editor-fold defaultstate="collapsed" desc="Private static final class constants">
    /**
     * Root logger
     */
    private static final Logger ROOT_LOGGER = LogManager.getLogManager().getLogger("");

    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(LogUtil.class.getName());
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Private static class variables">
    /**
     * File handler added to the root logger
     */
    private static FileHandler fileHandler;

    /**
     * Log file that is currently logged in
     */
    private static File logFile;
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Public static synchronized methods">
    /**
     * Configures the root logger called by all other loggers in the application not to print to console but to write
     * all logs to the log file specified in preferences. Also logs session start and sets as
     * default uncaught exception handler for all threads an object calling the root logger upon invocation. This method
     * should be invoked once upon starting the application. Then logging can be done by the individual class loggers
     * that will pass their logging messages to the root logger as default.
     * <p>
     * NOTE: Log file related methods need to be synchronized.
     *
     * @return true, if all actions were successfully executed; false, if an exception was thrown; in latter case the
     * root logger will be reset to default configuration
     */
    public static synchronized boolean initializeLoggingEnvironment() {
        try {
            //Configure root logger
            //Remove default console handler or a closed file handler after reset
            Handler[] tmpHandlers = LogUtil.ROOT_LOGGER.getHandlers();
            for (Handler tmpHandler : tmpHandlers) {
                LogUtil.ROOT_LOGGER.removeHandler(tmpHandler);
            }
            //Messages of levels INFO, WARNING and SEVERE will be logged only
            LogUtil.ROOT_LOGGER.setLevel(Level.INFO);
            String tmpLoggingDirectoryPathName = FileUtil.getAppDirPath() + File.separator
                    + BasicDefinitions.LOG_FILES_DIRECTORY + File.separator;
            File tmpLoggingDirectoryFile = new File(tmpLoggingDirectoryPathName);
            //If the directories do not exist already they are created
            if (!tmpLoggingDirectoryFile.exists()) {
                FileUtil.createDirectory(tmpLoggingDirectoryFile.getAbsolutePath());
            }
            String tmpLogFilePathName = tmpLoggingDirectoryPathName + BasicDefinitions.LOG_FILE_NAME
                    + "_"
                    + FileUtil.getTimeStampFileNameExtension();
            String tmpFinalLogFilePathName = FileUtil.getNonExistingFilePath(tmpLogFilePathName, BasicDefinitions.LOG_FILE_NAME_EXTENSION);
            File tmpLogFile = new File(tmpFinalLogFilePathName);
            boolean tmpFileWasCreated = FileUtil.createEmptyFile(tmpLogFile.getAbsolutePath());
            if (!tmpFileWasCreated) {
                throw new Exception("Log file " + tmpFinalLogFilePathName + " could not be created.");
            }
            if (!tmpLogFile.isFile() || !tmpLogFile.canWrite()) {
                throw new Exception("The designated log file " + tmpFinalLogFilePathName + " is not a file or can not be written to.");
            }
            LogUtil.logFile = tmpLogFile;
            LogUtil.fileHandler = new FileHandler(tmpFinalLogFilePathName, true);
            LogUtil.fileHandler.setFormatter(new SimpleFormatter());
            LogUtil.ROOT_LOGGER.addHandler(LogUtil.fileHandler);
            Thread.setDefaultUncaughtExceptionHandler(new Thread.UncaughtExceptionHandler() {
                @Override
                public void uncaughtException(Thread aThread, Throwable aThrowable) {
                    Logger.getLogger(aThread.getClass().getName()).log(Level.SEVERE, aThrowable.toString(), aThrowable);
                    if (aThread.getThreadGroup().getName().equals("main")) {
                        GuiUtil.GuiMessageAlert(Alert.AlertType.ERROR,
                                Message.get("Error.Notification.Title"),
                                null,
                                Message.get("Error.UnknownError"));
                        System.exit(-1);
                    }
                }
            });
            return true;
        } catch (Exception anException) {
            LogManager.getLogManager().reset();
            LogUtil.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            return false;
        }
    }

    /**
     * Reset log file. LogUtils.initializeLoggingEnvironment() will be called to reset the logging environment.
     * <p>
     * NOTE: Log file related methods need to be synchronized.
     * @return true, if all actions were successfully executed; false, if an exception was thrown or the log file in use
     * is not a file
     */
    public static synchronized boolean resetLogFile() {
        try {
            if (Objects.isNull(LogUtil.logFile) || !LogUtil.logFile.isFile()) {
                return false;
            }
            if (Objects.isNull(LogUtil.fileHandler)) {
                return false;
            }
            LogUtil.fileHandler.flush();
            LogUtil.fileHandler.close();
            boolean tmpFileWasDeleted = FileUtil.deleteSingleFile(LogUtil.logFile.getAbsolutePath());
            if (tmpFileWasDeleted) {
                boolean tmpWasLogEnvInitialized = LogUtil.initializeLoggingEnvironment();
                if (tmpWasLogEnvInitialized) {
                    LogUtil.LOGGER.info("Log file was reset.");
                }
                return true;
            } else {
                return false;
            }
        } catch (Exception anException) {
            LogUtil.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            return false;
        }
    }

    /**
     * Manages the folder the log-file get saved in if the folder exists.
     * If the folder holds more .txt files than a specific limit or a minimum of .txt files while exceeding
     * a maximum limit of bytes used, half of the .txt files get deleted and the method is called again.
     * Remaining LCK-files (suffix "*.txt.lck") are generally deleted out of the log-files' folder.
     * @author Samuel Behr
     */
    public static void manageLogFilesFolderIfExists() {
        Path tmpLogFileDirectory = Paths.get(FileUtil.getAppDirPath() + File.separator + BasicDefinitions.LOG_FILES_DIRECTORY);
        if (!(Files.exists(tmpLogFileDirectory) && Files.isDirectory(tmpLogFileDirectory))) {
            return;
        }
        //deleting all of the *.txt.lck files out of the log-files' folder
        try (DirectoryStream<Path> tmpLCKFilePaths = Files.newDirectoryStream(tmpLogFileDirectory, "*.txt.lck")) {
            for (Path tmpLCKFilePath : tmpLCKFilePaths) {
                try {
                    Files.delete(tmpLogFileDirectory.resolve(tmpLCKFilePath));
                } catch (IOException anException) {
                    //TODO
                }
            }
        } catch (IOException anException) {
            //TODO
        }
        File [] tmpLogFiles = tmpLogFileDirectory.toFile().listFiles((dir, name) -> name.endsWith(".txt"));
        int tmpTotalOfBytesUsed = 0;
        for (File tmpLogFile : tmpLogFiles) {
            tmpTotalOfBytesUsed += tmpLogFile.length();
        }
        //managing the log-files if the limits are exceeded
        //the parameters of this if statement's condition should be changed with caution or otherwise an infinite loop is risked
        if (tmpLogFiles.length > BasicDefinitions.UPPER_LIMIT_OF_LOG_FILES || (tmpTotalOfBytesUsed > BasicDefinitions.LIMIT_OF_BYTES_USED_BY_LOG_FILES
                && tmpLogFiles.length > BasicDefinitions.LOWER_LIMIT_OF_LOG_FILES)) {
            Arrays.sort(tmpLogFiles, Comparator.comparingLong(File::lastModified));
            //deleting the first half of the files of the sorted File array out of the log-files' folder
            int tmpDeletedLogFilesCounter = 0;  //TODO: remove??
            for (int i = 0; i < (tmpLogFiles.length * BasicDefinitions.FACTOR_TO_TRIM_LOG_FILE_FOLDER); i++) {
                try {
                    Files.delete(tmpLogFileDirectory.resolve(tmpLogFiles[i].toPath()));
                    tmpDeletedLogFilesCounter++;
                } catch (IOException anException) {
                    //TODO
                }
            }
            System.out.println("Deleted files counter: " + tmpDeletedLogFilesCounter);  //TODO: remove!
            //calling the method again to check whether the limits are still exceeded
            LogUtil.manageLogFilesFolderIfExists();
            //log-files' folder has been managed and ... log-files have been deleted. //TODO??
        }
    }
    // </editor-fold>
}
