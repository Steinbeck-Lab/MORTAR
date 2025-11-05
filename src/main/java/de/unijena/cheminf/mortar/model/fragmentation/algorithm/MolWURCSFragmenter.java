/*
 * MORTAR - MOlecule fRagmenTAtion fRamework
 * Copyright (C) 2025  Felix Baensch, Jonas Schaub (felix.j.baensch@gmail.com, jonas.schaub@uni-jena.de)
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

package de.unijena.cheminf.mortar.model.fragmentation.algorithm;

import de.unijena.cheminf.mortar.message.Message;
import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.util.BasicDefinitions;
import de.unijena.cheminf.mortar.model.util.CollectionUtil;

import javafx.beans.property.Property;
import javafx.beans.property.SimpleBooleanProperty;

import org.glycoinfo.MolWURCS.exchange.fromWURCS.WURCSGraphToMolecule;
import org.glycoinfo.MolWURCS.io.WURCSWriter;
import org.glycoinfo.MolWURCS.util.analysis.HighEnergySiteFinder;
import org.glycoinfo.WURCSFramework.util.WURCSException;
import org.glycoinfo.WURCSFramework.util.WURCSFactory;
import org.glycoinfo.WURCSFramework.wurcs.graph.WURCSGraph;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Wrapper class that makes sugar extraction via the
 * <a href="https://pubs.acs.org/doi/full/10.1021/ci400571e">WURCS glycoside identifier</a>
 * available in MORTAR, using the <a href="https://doi.org/10.1007/s00216-024-05508-1">MolWURCS</a>
 * implementation of the identifier.
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class MolWURCSFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm used in this fragmenter.
     */
    public static final String ALGORITHM_NAME = "WURCS algorithm";

    /**
     * Default value for whether the aglycone part should be included in the output fragments.
     */
    public static final boolean OUTPUT_WITH_AGLYCONE_SETTING_DEFAULT = false;

    /**
     * Default value for whether double check should be performed in WURCSWriter.
     */
    public static final boolean DOUBLE_CHECK_SETTING_DEFAULT = false;

    /**
     * Default value for whether WURCS normalization should be applied before retranslating WURCS to molecule.
     */
    public static final boolean NORMALIZE_WURCS_BEFORE_RETRANSLATION_DEFAULT = true;
    //</editor-fold>
    //
    //<editor-fold desc="Private final variables">
    /**
     * All settings of this fragmenter, encapsulated in JavaFX properties for binding in GUI.
     */
    private final List<Property<?>> settings;

    /**
     * Map to store pairs of {@literal <setting name, tooltip text>}.
     */
    private final HashMap<String, String> settingNameTooltipTextMap;

    /**
     * Map to store pairs of {@literal <setting name, display name>}.
     */
    private final HashMap<String, String> settingNameDisplayNameMap;

    /**
     * Property for whether the aglycone part should be included in the output fragments.
     */
    private final SimpleBooleanProperty outputWithAglyconeSetting;

    /**
     * Property for whether double check should be performed in WURCSWriter.
     */
    private final SimpleBooleanProperty doubleCheckSetting;

    /**
     * Property for whether WURCS normalization should be applied before retranslating WURCS to molecule.
     */
    private final SimpleBooleanProperty normalizeWURCSBeforeRetranslationSetting;
    //</editor-fold>
    //
    //<editor-fold desc="Private static final constants">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(MolWURCSFragmenter.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor, all settings are initialised with their default values.
     */
    public MolWURCSFragmenter() {
        int tmpNumberOfSettings = 3;
        this.settings = new ArrayList<>(tmpNumberOfSettings);
        int tmpInitialCapacityForSettingMaps = CollectionUtil.calculateInitialHashCollectionCapacity(
                tmpNumberOfSettings,
                BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameTooltipTextMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.settingNameDisplayNameMap = new HashMap<>(tmpInitialCapacityForSettingMaps, BasicDefinitions.DEFAULT_HASH_COLLECTION_LOAD_FACTOR);
        this.outputWithAglyconeSetting = new SimpleBooleanProperty(this, "Output with aglycone setting",
                MolWURCSFragmenter.OUTPUT_WITH_AGLYCONE_SETTING_DEFAULT);
        this.settings.add(this.outputWithAglyconeSetting);
        this.settingNameDisplayNameMap.put(this.outputWithAglyconeSetting.getName(),
                Message.get("MolWURCSFragmenter.outputWithAglyconeSetting.displayName"));
        this.settingNameTooltipTextMap.put(this.outputWithAglyconeSetting.getName(),
                Message.get("MolWURCSFragmenter.outputWithAglyconeSetting.tooltip"));
        this.doubleCheckSetting = new SimpleBooleanProperty(this, "Double check setting",
                MolWURCSFragmenter.DOUBLE_CHECK_SETTING_DEFAULT);
        this.settings.add(this.doubleCheckSetting);
        this.settingNameDisplayNameMap.put(this.doubleCheckSetting.getName(),
                Message.get("MolWURCSFragmenter.doubleCheckSetting.displayName"));
        this.settingNameTooltipTextMap.put(this.doubleCheckSetting.getName(),
                Message.get("MolWURCSFragmenter.doubleCheckSetting.tooltip"));
        this.normalizeWURCSBeforeRetranslationSetting = new SimpleBooleanProperty(this,
                "Normalize WURCS before retranslation setting",
                MolWURCSFragmenter.NORMALIZE_WURCS_BEFORE_RETRANSLATION_DEFAULT);
        this.settings.add(this.normalizeWURCSBeforeRetranslationSetting);
        this.settingNameDisplayNameMap.put(this.normalizeWURCSBeforeRetranslationSetting.getName(),
                Message.get("MolWURCSFragmenter.normalizeWURCSBeforeRetranslationSetting.displayName"));
        this.settingNameTooltipTextMap.put(this.normalizeWURCSBeforeRetranslationSetting.getName(),
                Message.get("MolWURCSFragmenter.normalizeWURCSBeforeRetranslationSetting.tooltip"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties get">
    /**
     * Returns the currently set value of the output with aglycone setting.
     *
     * @return true, if the fragments should contain the aglycone part as well; false, if only the sugar moiety should be output
     */
    public boolean getOutputWithAglyconeSetting() {
        return this.outputWithAglyconeSetting.get();
    }

    /**
     * Returns the currently set value of the double check setting.
     *
     * @return true, if double check is enabled; false, otherwise
     */
    public boolean getDoubleCheckSetting() {
        return this.doubleCheckSetting.get();
    }

    /**
     * Returns the currently set value of the 'normalize WURCS before retranslation setting'.
     *
     * @return true, if WURCS normalization is applied before retranslating WURCS to molecule; false, otherwise
     */
    public boolean getNormalizeWURCSBeforeRetranslationSetting() {
        return this.normalizeWURCSBeforeRetranslationSetting.get();
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public properties set">
    /**
     * Sets the value of the output with aglycone setting.
     *
     * @param aBoolean true, if the fragments should contain the aglycone part as well; false, if only the sugar moiety should be output
     */
    public void setOutputWithAglyconeSetting(boolean aBoolean) {
        this.outputWithAglyconeSetting.set(aBoolean);
    }

    /**
     * Sets the value of the double check setting.
     *
     * @param aBoolean true, if double check should be performed by the WURCSWriter; false, otherwise
     */
    public void setDoubleCheckSetting(boolean aBoolean) {
        this.doubleCheckSetting.set(aBoolean);
    }

    /**
     * Sets the value of the 'normalize WURCS before retranslation' setting.
     *
     * @param aBoolean true, if WURCS normalization should be applied before retranslating WURCS to molecule; false, otherwise
     */
    public void setNormalizeWURCSBeforeRetranslationSetting(boolean aBoolean) {
        this.normalizeWURCSBeforeRetranslationSetting.set(aBoolean);
    }
    //</editor-fold>
    //
    //<editor-fold desc="IMoleculeFragmenter methods">
    //without the empty line, the code folding does not work properly here...

    @Override
    public List<Property<?>> settingsProperties() {
        return this.settings;
    }

    @Override
    public Map<String, String> getSettingNameToTooltipTextMap() {
        return this.settingNameTooltipTextMap;
    }

    @Override
    public Map<String, String> getSettingNameToDisplayNameMap() {
        return this.settingNameDisplayNameMap;
    }

    @Override
    public String getFragmentationAlgorithmName() {
        return MolWURCSFragmenter.ALGORITHM_NAME;
    }

    @Override
    public String getFragmentationAlgorithmDisplayName() {
        return Message.get("MolWURCSFragmenter.displayName");
    }

    @Override
    public IMoleculeFragmenter copy() {
        MolWURCSFragmenter tmpCopy = new MolWURCSFragmenter();
        tmpCopy.setOutputWithAglyconeSetting(this.getOutputWithAglyconeSetting());
        tmpCopy.setDoubleCheckSetting(this.getDoubleCheckSetting());
        tmpCopy.setNormalizeWURCSBeforeRetranslationSetting(this.getNormalizeWURCSBeforeRetranslationSetting());
        return tmpCopy;
    }

    @Override
    public void restoreDefaultSettings() {
        this.outputWithAglyconeSetting.set(MolWURCSFragmenter.OUTPUT_WITH_AGLYCONE_SETTING_DEFAULT);
        this.doubleCheckSetting.set(MolWURCSFragmenter.DOUBLE_CHECK_SETTING_DEFAULT);
        this.normalizeWURCSBeforeRetranslationSetting.set(MolWURCSFragmenter.NORMALIZE_WURCS_BEFORE_RETRANSLATION_DEFAULT);
    }

    @Override
    public List<IAtomContainer> fragmentMolecule(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Parameter tests">
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpCanBeFragmented = this.canBeFragmented(aMolecule);
        if (!tmpCanBeFragmented) {
            throw new IllegalArgumentException("Given molecule cannot be fragmented but should be filtered or preprocessed first.");
        }
        //</editor-fold>
        //MolWURCS applies preprocessing to atom container, so safer to work with a clone here
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        //output writer for WURCSWriter
        StringWriter tmpStringWriter = new StringWriter();
        //note: safer to always initialise a new writer
        WURCSWriter tmpWURCSWriter = new WURCSWriter(tmpStringWriter);
        tmpWURCSWriter.setOutputWithAglycone(this.outputWithAglyconeSetting.get());
        tmpWURCSWriter.setDoDoubleCheck(this.doubleCheckSetting.get());
        //set property key for molecule title for log messages
        tmpWURCSWriter.setTitlePropertyID(Importer.MOLECULE_NAME_PROPERTY_KEY);
        //little preprocessing to prevent WURCS from trying to deduce stereo config from coordinates that are not there
        for (IBond bond : tmpMoleculeClone.bonds()) {
            if (bond.getDisplay() == IBond.Display.Solid && bond.getOrder() == IBond.Order.DOUBLE) {
                bond.setDisplay(IBond.Display.Crossed);
            }
        }
        //preprocessing necessary for aromaticity perception in WURCS preprocessing
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
        } catch (CDKException anException) {
            throw new IllegalArgumentException("Could not perceive atom types and configure atoms of the molecule to fragment("
                    + tmpMoleculeClone.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY) + ").", anException);
        }
        //writes the WURCS representation to the StringWriter
        tmpWURCSWriter.writeAtomContainer(tmpMoleculeClone);
        try {
            tmpWURCSWriter.close();
        } catch (IOException anException) {
            MolWURCSFragmenter.LOGGER.log(Level.WARNING, String.format("Could not close WURCSWriter after writing molecule: %s",
                    tmpMoleculeClone.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)), anException);
        }
        //retrieve generated WURCS String
        String tmpWURCSString = tmpStringWriter.toString();
        //parse WURCS String back to Molecule(s)
        // if multiple sugar moieties were detected, they are separated by new lines
        List<String> tmpSeparateSugarWURCSCodesList = tmpWURCSString.lines().toList();
        if (tmpSeparateSugarWURCSCodesList.isEmpty() || (tmpSeparateSugarWURCSCodesList.size() == 1
                && (Objects.isNull(tmpSeparateSugarWURCSCodesList.getFirst()) || tmpSeparateSugarWURCSCodesList.getFirst().isBlank()))) {
            return new ArrayList<>(0);
        }
        List<IAtomContainer> tmpFragments = new ArrayList<>(tmpSeparateSugarWURCSCodesList.size());
        WURCSGraphToMolecule tmpWURCSGraphToMol = new WURCSGraphToMolecule();
        for (String tmpWURCSCode : tmpSeparateSugarWURCSCodesList) {
            //if there was an issue, the WURCSWriter might have logged the exception to console and returned an empty line
            // we want to detect such cases and properly throw an exception, so that it gets counted for the log
            if (Objects.isNull(tmpWURCSCode) || tmpWURCSCode.isBlank()) {
                throw new NullPointerException("Error while creating WURCS code: " + tmpWURCSCode
                        + " molecule name: " + tmpMoleculeClone.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY));
            }
            WURCSFactory tmpWURCSFactory;
            try {
                //factory has to be newly instantiated for every molecule
                tmpWURCSFactory = new WURCSFactory(tmpWURCSCode, this.normalizeWURCSBeforeRetranslationSetting.get());
                WURCSGraph tmpWURCSGraph = tmpWURCSFactory.getGraph();
                try {
                    tmpWURCSGraphToMol.start(tmpWURCSGraph);
                } catch (NullPointerException anException) {
                    //re-throwing to add the molecule name to the message
                    throw new NullPointerException("Error while parsing WURCS code: " + tmpWURCSCode
                            + " molecule name: " + tmpMoleculeClone.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY)
                            + " original message: " + anException.getMessage());
                }
                IAtomContainer tmpMolecule = tmpWURCSGraphToMol.getMolecule();
                if (!Objects.isNull(tmpMolecule)) {
                    tmpFragments.add(tmpMolecule);
                } else {
                    throw new NullPointerException("Could not retranslate WURCS code to molecule: " + tmpWURCSCode
                            + " molecule name: " + tmpMoleculeClone.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY));
                }
            } catch (WURCSException tmpException) {
                throw new IllegalArgumentException("Error while parsing WURCS code: " + tmpWURCSCode
                        + " molecule name: " + tmpMoleculeClone.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY),
                        tmpException);
                //continue;
            }
        }
        return tmpFragments;
    }

    @Override
    public boolean shouldBeFiltered(IAtomContainer aMolecule) {
        if (Objects.isNull(aMolecule) || aMolecule.isEmpty()) {
            return true;
        }
        //note: we could check for many other things that WURCS does not support here, but it checks again on its own anyway...
        //checks for radicals, high-energy bonds and rings
        return HighEnergySiteFinder.find(aMolecule);
    }

    @Override
    public boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        //preprocessing is done internally by WURCS later
        return false;
    }

    @Override
    public boolean canBeFragmented(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered = this.shouldBeFiltered(aMolecule);
        boolean tmpShouldBePreprocessed = this.shouldBePreprocessed(aMolecule);
        return !(tmpShouldBeFiltered || tmpShouldBePreprocessed);
    }

    @Override
    public IAtomContainer applyPreprocessing(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered = this.shouldBeFiltered(aMolecule);
        if (tmpShouldBeFiltered) {
            throw new IllegalArgumentException("The given molecule cannot be preprocessed but should be filtered.");
        }
        //preprocessing is done internally by WURCS later
        return aMolecule.clone();
    }
}
