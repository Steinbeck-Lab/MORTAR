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

import javafx.beans.property.Property;

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

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * TODO
 */
public class MolWURCSFragmenter implements IMoleculeFragmenter {
    //<editor-fold desc="Public static final constants">
    /**
     * Name of the algorithm used in this fragmenter.
     */
    public static final String ALGORITHM_NAME = "WURCS algorithm";
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
        this.settings = List.of();
        this.settingNameTooltipTextMap = new HashMap<>();
        this.settingNameDisplayNameMap = new HashMap<>();
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
        return tmpCopy;
    }

    @Override
    public void restoreDefaultSettings() {

    }

    @Override
    public List<IAtomContainer> fragmentMolecule(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        StringWriter tmpStringWriter = new StringWriter();
        WURCSWriter tmpWURCSWriter = new WURCSWriter(tmpStringWriter);
        tmpWURCSWriter.setOutputWithAglycone(false); //todo this can be a setting
        tmpWURCSWriter.setDoDoubleCheck(true); //todo true? - setting?
        tmpWURCSWriter.setTitlePropertyID(Importer.MOLECULE_NAME_PROPERTY_KEY);
        //little preprocessing to prevent WURCS from trying to deduce stereo config from coordinates that are not there
        for (IBond bond : aMolecule.bonds()) {
            if (bond.getDisplay() == IBond.Display.Solid && bond.getOrder() == IBond.Order.DOUBLE) {
                bond.setDisplay(IBond.Display.Crossed);
            }
        }
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        } catch (CDKException anException) {
            throw new IllegalArgumentException("Could not perceive atom types and configure atoms of the molecule to fragment("
                    + aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY) + ").", anException);
        }
        //writes the WURCS representation to the StringWriter
        tmpWURCSWriter.writeAtomContainer(aMolecule);
        // Get WURCS String
        String tmpWURCSString = tmpStringWriter.toString();
        // Parse WURCS String back to Molecule(s)
        // if multiple sugar moieties were detected, they are separated by new lines
        List<String> tmpSeparateSugarWURCSCodesList = tmpWURCSString.lines().toList();
        List<IAtomContainer> tmpFragments = new ArrayList<>(tmpSeparateSugarWURCSCodesList.size());
        for (String tmpWURCSCode : tmpSeparateSugarWURCSCodesList) {
            WURCSFactory tmpWURCSFactory = null;
            try {
                tmpWURCSFactory = new WURCSFactory(tmpWURCSCode, true); //todo investigate param, setting?
            } catch (WURCSException tmpException) {
                MolWURCSFragmenter.LOGGER.log(Level.SEVERE, "Error while parsing WURCS code: " + tmpWURCSCode + " molecule name: " + aMolecule.getProperty(Importer.MOLECULE_NAME_PROPERTY_KEY), tmpException);
                continue;
            }
            WURCSGraph tmpWURCSGraph = tmpWURCSFactory.getGraph();
            WURCSGraphToMolecule tmpWURCSGraphToMol = new WURCSGraphToMolecule();
            tmpWURCSGraphToMol.start(tmpWURCSGraph);
            tmpFragments.add(tmpWURCSGraphToMol.getMolecule());
        }
        return tmpFragments;
    }

    @Override
    public boolean shouldBeFiltered(IAtomContainer aMolecule) {
        if (aMolecule.isEmpty()) {
            return true;
        }
        //note: we could check for many other things that WURCS does not support here but it checks again on its own again anyway...
        return HighEnergySiteFinder.find(aMolecule);
    }

    @Override
    public boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        //done internally by WURCS
        return false;
    }

    @Override
    public boolean canBeFragmented(IAtomContainer aMolecule) throws NullPointerException {
        return !this.shouldBeFiltered(aMolecule) && this.shouldBePreprocessed(aMolecule);
    }

    @Override
    public IAtomContainer applyPreprocessing(IAtomContainer aMolecule) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //do nothing
        return aMolecule;
    }
}
