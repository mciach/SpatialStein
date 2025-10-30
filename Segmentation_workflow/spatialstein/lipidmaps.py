def load_lipidmaps_formulas(lipidmaps_file_path):
    """
    A helper function to load lipidmaps formulas from a tab-separated file.
    We recommend preparing a custom file with lipids from selected classes
    or using the file provided with spatialstein.  
    Selects unique formulas with elements CHNOP.
    Returns a list of molecular formulas and the mapping to their
    IDs, classes and subclasses.
    Note that a single molecular formula may correspond to multiple lipids,
    in which case each of their class is returned.   
    """
    lipid_formulas = []
    formula_to_subclasses = {}
    formula_to_classes = {}
    formula_to_IDs = {}
    with open(lipidmaps_file_path) as h:
        for l in h:
            l = l.strip().split('\t')
            formula = l[3] 
            lipid_class = l[1]
            lipid_subclass = l[2]
            lipid_ID = l[0]
            if set(formula).issubset(set('CHNOP0123456789')): # skipping atypical formulas, e.g. with Deuterium
                lipid_formulas.append(formula)
                try:
                    formula_to_subclasses[formula].append(lipid_subclass)
                except KeyError:
                    formula_to_subclasses[formula] = [lipid_subclass]
                try:
                    formula_to_classes[formula].append(lipid_class)
                except KeyError:
                    formula_to_classes[formula] = [lipid_class]
                try:
                    formula_to_IDs[formula].append(lipid_ID)
                except KeyError:
                    formula_to_IDs[formula] = [lipid_ID]
    lipid_formulas = list(set(lipid_formulas))
    nb_of_classes = len(set([x for formula in formula_to_classes for x in formula_to_classes[formula]]))
    nb_of_subclasses = len(set([x for formula in formula_to_subclasses for x in formula_to_subclasses[formula]]))
    nb_of_IDs = len(set([x for formula in formula_to_IDs for x in formula_to_IDs[formula]]))
    print('Loaded %i unique lipid formulas' % len(lipid_formulas))
    print('Loaded %i unique lipid classes' % nb_of_classes)
    print('Loaded %i unique lipid subclasses' % nb_of_subclasses)
    print('Loaded %i unique lipid IDs' % nb_of_IDs)
    return(lipid_formulas, formula_to_IDs, formula_to_subclasses, formula_to_classes)
