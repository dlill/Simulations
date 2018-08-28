# Description of the subfolders/models



* Cascade
    * The standard three tier cascade with gainuf/gainlf feedbacks
    * synthesis and degradation
    * Used for Figure 4

* Cascade1long1short
    * The Cascade with (too) many Feedback-complexes, some of the reactions don't make sense. See the Figure "Copy of cascade1long1short.xml" where the same complexes are in color. They can only be in one pathway, else it's nonsensical
    * The too many complexes make the whole thing really not nice to handle

* CascadeBondGraph
    * Attempt to reduce the complexity of Cascade1long1short

* CascadeComprehensiveFeedbacks
    * Cascade with one long and two short feedbacks which are all modeled via complexes
    * Analysis file (.Rmd) is in the folder itself
        * Analysis happens via Shiny gadget
    * Analytical structural investigation failed, but numerically, it works, they are just the numerical sensitivities (i.e. Global responses)
        * The global responses have different signs for different combinations, therefore it is clear that with the free parameters, it is possibly to minimize the different connections

* Directcomplexfeedback
    *

* Hub
    *

* HubInBetween
    *

* HubInCascade
    *

* Inactivebinding
    * Inactive forms binding without modifying each other. 
    * Method doesn't work here

* MassActionCascade2Modules
    *

* MST2Raf
    *

* MST2Raf all the other modules
    *

* Phosphatase
    *

* Prabakaran
    *

* Prabakaranphosphatase
    *

* RafMekErk
    *

* Ras
    *

* RegulatoryCascade
    *

* ScaffoldCascade
    *
