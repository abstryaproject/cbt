/*************************************************
 * 1️⃣ ADD YOUR QUESTIONS HERE (MANUALLY)
 * ➤ Each question must have:
 *    - question: string
 *    - options: array of strings
 *    - answer: index of correct option (0-based)
 *    - explanation: string
 *    - subject: string (e.g., "Math", "Physics", "English")
 *************************************************/

let quizData1 = [
  {
    question: "Which molecular marker is most commonly used for bacterial phylogenetic analysis?",
    options: ["16S rRNA gene", "18S rRNA gene", "ITS region", "COI gene"],
    answer: 0,
    explanation: "16S rRNA is highly conserved in bacteria and widely used for phylogenetic studies."
  },
  {
    question: "Which staining method differentiates Gram-positive from Gram-negative bacteria?",
    options: ["Acid-fast stain", "Gram stain", "Endospore stain", "Capsule stain"],
    answer: 1,
    explanation: "Gram stain distinguishes bacteria based on cell wall peptidoglycan thickness."
  },
  {
    question: "Which genus is known for forming endospores and belongs to Firmicutes?",
    options: ["Bacillus", "Escherichia", "Streptomyces", "Pseudomonas"],
    answer: 0,
    explanation: "Bacillus species produce endospores and are members of Firmicutes."
  },
  {
    question: "Which molecular technique allows identification of unculturable microorganisms directly from environmental samples?",
    options: ["Culture on selective media", "16S rRNA gene sequencing", "Gram staining", "Electron microscopy"],
    answer: 1,
    explanation: "Sequencing 16S rRNA directly from environmental DNA detects unculturable microbes."
  },
  {
    question: "Which group of archaea is commonly associated with extreme halophilic environments?",
    options: ["Euryarchaeota", "Crenarchaeota", "Nanoarchaeota", "Korarchaeota"],
    answer: 0,
    explanation: "Halophilic archaea are mainly Euryarchaeota adapted to high salt concentrations."
  },
  {
    question: "Which taxonomic rank is immediately above genus in the hierarchical classification?",
    options: ["Family", "Order", "Class", "Species"],
    answer: 0,
    explanation: "The order is above family, family is above genus."
  },
  {
    question: "Which feature distinguishes fungi from bacteria?",
    options: ["Presence of peptidoglycan", "Presence of chitin in cell wall", "Nucleoid DNA", "Capsule formation"],
    answer: 1,
    explanation: "Fungal cell walls contain chitin, whereas bacteria have peptidoglycan."
  },
  {
    question: "Which virus classification system is based on nucleic acid type, strandness, and replication strategy?",
    options: ["Linnaean system", "Baltimore classification", "ICTV system", "Gram classification"],
    answer: 1,
    explanation: "Baltimore classification groups viruses by genome type and replication."
  },
  {
    question: "Which genus is a Gram-negative, photosynthetic, and nitrogen-fixing bacterium?",
    options: ["Rhizobium", "Rhodobacter", "Azotobacter", "Clostridium"],
    answer: 1,
    explanation: "Rhodobacter is Gram-negative, photosynthetic, and capable of nitrogen fixation."
  },
  {
    question: "Which molecular method is preferred for high-resolution identification at the species level among fungi?",
    options: ["16S rRNA sequencing", "ITS region sequencing", "Gram staining", "Lipid profiling"],
    answer: 1,
    explanation: "The Internal Transcribed Spacer (ITS) region is highly variable and used for fungal species identification."
  },
  {
    question: "Which bacterial phylum contains members that are extremophiles, such as thermophilic sulfur reducers?",
    options: ["Proteobacteria", "Firmicutes", "Crenarchaeota", "Actinobacteria"],
    answer: 2,
    explanation: "Crenarchaeota includes thermophilic archaea adapted to high temperature and sulfur-rich environments."
  },
  {
    question: "Which of the following is an example of a polyphasic approach in microbial taxonomy?",
    options: [
      "Using only 16S rRNA sequencing",
      "Combining morphology, physiology, and molecular data",
      "Culture on selective media only",
      "Observation under light microscope"
    ],
    answer: 1,
    explanation: "Polyphasic taxonomy integrates phenotypic, chemotaxonomic, and genotypic data for classification."
  },
  {
    question: "Which genus of actinobacteria is well-known for producing antibiotics?",
    options: ["Streptomyces", "Bacillus", "Clostridium", "Pseudomonas"],
    answer: 0,
    explanation: "Streptomyces species produce many clinically important antibiotics."
  },
  {
    question: "Which of the following is a limitation of phenotypic methods in microbial taxonomy?",
    options: ["Cannot detect unculturable microbes", "Too fast to interpret", "Highly expensive", "Requires PCR"],
    answer: 0,
    explanation: "Phenotypic methods rely on culturable microbes and fail to detect unculturable organisms."
  },
  {
    question: "Which viral family includes double-stranded RNA viruses?",
    options: ["Reoviridae", "Retroviridae", "Coronaviridae", "Adenoviridae"],
    answer: 0,
    explanation: "Reoviridae contains dsRNA genomes."
  },
  {
    question: "Which taxonomic tool uses DNA-DNA hybridization to determine genomic similarity?",
    options: ["16S rRNA sequencing", "Whole genome sequencing", "DDH", "PCR-RFLP"],
    answer: 2,
    explanation: "DNA-DNA hybridization quantifies genome-level similarity between strains."
  },
  {
    question: "Which genus of cyanobacteria is filamentous and heterocyst-forming?",
    options: ["Nostoc", "Synechococcus", "Prochlorococcus", "Anabaena"],
    answer: 3,
    explanation: "Anabaena forms heterocysts for nitrogen fixation."
  },
  {
    question: "Which type of microbial taxonomy is based on evolutionary relationships?",
    options: ["Phenetic", "Phylogenetic", "Ecological", "Morphological"],
    answer: 1,
    explanation: "Phylogenetic classification reflects evolutionary lineage."
  },
  {
    question: "Which technique can rapidly classify microbes without culturing by analyzing rRNA sequences?",
    options: ["Metagenomics", "Culture plating", "Gram staining", "Endospore staining"],
    answer: 0,
    explanation: "Metagenomics sequences rRNA genes directly from environmental samples, detecting unculturable microbes."
  }
];

let quizData2 = [
  {
    question: "Which bacterium is most often genetically engineered for high‑yield production of heterologous proteins in industry?",
    options: ["Bacillus subtilis", "Escherichia coli", "Clostridium perfringens", "Pseudomonas aeruginosa"],
    answer: 1,
    explanation: "E. coli is widely used for recombinant protein production due to fast growth, defined genetics, and many vector tools."
  },
  {
    question: "What is the primary role of *Saccharomyces cerevisiae* in industrial fermentation?",
    options: ["Lactic acid production", "Cellulase secretion", "Ethanol production", "Antibiotic biosynthesis"],
    answer: 2,
    explanation: "S. cerevisiae ferments sugars anaerobically to ethanol, a key process in brewing and biofuel production."
  },
  {
    question: "In the activated sludge process, the term 'food‑to‑microorganism ratio (F/M)' influences:",
    options: ["Biomass settling", "Sludge age", "Oxygen uptake", "All of the above"],
    answer: 3,
    explanation: "F/M ratio affects microbial growth, settling characteristics, oxygen demand, and treatment efficiency."
  },
  {
    question: "Which group of microbes is primarily responsible for nitrification?",
    options: ["Denitrifiers", "Nitrogen fixers", "Ammonia oxidizers", "Methanogens"],
    answer: 2,
    explanation: "Nitrosomonas and Nitrobacter catalyze ammonia → nitrite → nitrate in aerobic conditions."
  },
  {
    question: "Which compound is a volatile fatty acid important in anaerobic digestion?",
    options: ["Acetate", "Glycerol", "Urea", "Glucose"],
    answer: 0,
    explanation: "Acetate is a key intermediate in anaerobic digestion, feeding methanogens for methane production."
  },
  {
    question: "What is the primary advantage of using *Aspergillus niger* over bacterial producers for citric acid?",
    options: ["Higher growth rate", "Lower pH tolerance", "Spore‑free fermentation", "No need oxygen"],
    answer: 1,
    explanation: "A. niger can tolerate very low pH, reducing contamination and improving citric acid yields."
  },
  {
    question: "In microbial fuel cells, electrons originate from:",
    options: ["Photosynthesis", "Organic substrate oxidation", "Nitrogen fixation", "Methanogenesis"],
    answer: 1,
    explanation: "Microbes oxidize organics and transfer electrons to electrodes, generating current."
  },
  {
    question: "Which analytical technique is used to monitor microbial community shifts in wastewater systems?",
    options: ["PCR‑DGGE", "Gram staining", "Endospore count", "pH meter"],
    answer: 0,
    explanation: "PCR‑DGGE (Denaturing Gradient Gel Electrophoresis) profiles community changes at the DNA level."
  },
  {
    question: "Which of the following is a major limitation of composting?",
    options: ["High cost", "Slow degradation of lignin", "Methane emission", "No pathogen reduction"],
    answer: 1,
    explanation: "Lignin is recalcitrant and slows compost maturation; other processes handle it poorly."
  },
  {
    question: "Recombinant *Bacillus thuringiensis* is used for industrial production of:",
    options: ["Biofertilizers", "Biopesticides", "Bioethanol", "Amino acids"],
    answer: 1,
    explanation: "Bt produces insecticidal proteins; recombinant strains improve yield and specificity."
  },
  {
    question: "Which enzyme class is most important for starch processing in industry?",
    options: ["Lipases", "Amylases", "Proteases", "Cellulases"],
    answer: 1,
    explanation: "Amylases hydrolyze starch into sugars, essential in brewing and sweetener production."
  },
  {
    question: "During nitrification, the key electron donor is:",
    options: ["Ammonia", "Nitrate", "Nitrogen gas", "Organic carbon"],
    answer: 0,
    explanation: "Ammonia oxidation releases electrons used by nitrifiers."
  },
  {
    question: "Which of the following best describes bioaugmentation?",
    options: [
      "Adding enzymes to accelerate composting",
      "Introducing specialized microbes to contaminated sites",
      "Using chemicals to change pH",
      "Sterilizing wastewater"
    ],
    answer: 1,
    explanation: "Bioaugmentation adds microbes that metabolize specific pollutants to improve cleanup."
  },
  {
    question: "What is the microbial basis of cheese flavor differentiation?",
    options: ["Lactic acid alone", "Protein proteolysis by microbial enzymes", "Salt diffusion", "Temperature"],
    answer: 1,
    explanation: "Proteolysis by microbes and their enzymes releases peptides and amino acids that produce distinct flavors."
  },
  {
    question: "Which microbial group fixes nitrogen in aerobic soils?",
    options: ["Cyanobacteria", "Rhizobia", "Clostridia", "Thiobacilli"],
    answer: 1,
    explanation: "Rhizobia form symbioses with legumes and fix nitrogen even in aerobic soil microsites."
  },
  {
    question: "The rate‑limiting step in anaerobic digestion is:",
    options: ["Hydrolysis", "Acetogenesis", "Methanogenesis", "Denitrification"],
    answer: 0,
    explanation: "Hydrolysis of complex polymers into monomers is often the slowest step in anaerobic digestion."
  },
  {
    question: "Which compound indicates effective pathogen reduction in treated wastewater?",
    options: ["High COD", "Low BOD", "Coliform count reduction", "Increase in ammonia"],
    answer: 2,
    explanation: "Reduction in coliform bacteria indicates successful disinfection."
  },
  {
    question: "Which industrial microorganism produces penicillin?",
    options: ["Penicillium chrysogenum", "Penicillium roqueforti", "Aspergillus flavus", "Streptomyces griseus"],
    answer: 0,
    explanation: "P. chrysogenum is the principal commercial penicillin producer."
  },
  {
    question: "Metagenomics studies in environmental microbiology are primarily used to:",
    options: [
      "Count colony forming units",
      "Sequence DNA from mixed communities",
      "Increase biomass yield",
      "Measure pH"
    ],
    answer: 1,
    explanation: "Metagenomics sequences environmental DNA to characterize unculturable microbes."
  },
  {
    question: "Which is the most energy‑efficient stage of wastewater treatment for reducing organic load?",
    options: ["Primary sedimentation", "Activated sludge", "Tertiary filtration", "Disinfection"],
    answer: 1,
    explanation: "Activated sludge biologically oxidizes organic matter efficiently with high removal rates."
  },
  {
    question: "In industrial fermentations, foam formation is usually controlled by:",
    options: ["Antifoam agents", "Lower temperature", "Increased agitation", "UV irradiation"],
    answer: 0,
    explanation: "Antifoams reduce surface tension and prevent excessive foam."
  },
  {
    question: "Decomposition of lignocellulose in compost is primarily carried out by:",
    options: ["Fungi", "Protozoa", "Methanogens", "Viruses"],
    answer: 0,
    explanation: "Fungi secrete enzymes that degrade lignin and cellulose."
  },
  {
    question: "Which of these is a dissimilatory nitrate reduction process producing gaseous nitrogen?",
    options: ["Nitrogen fixation", "Assimilatory nitrate reduction", "Denitrification", "Nitritation"],
    answer: 2,
    explanation: "Denitrification reduces nitrates to N2 or N2O under anaerobic conditions."
  },
  {
    question: "Microbial spoilage in fermented foods is often prevented by:",
    options: ["High oxygen", "Low pH", "High temperature", "Low salt"],
    answer: 1,
    explanation: "Low pH inhibits many spoilage microbes while favoring desired fermenters."
  },
  {
    question: "Which fermentation product is key in sourdough bread flavor and texture?",
    options: ["Lactic acid", "Acetic acid", "Ethanol", "All of the above"],
    answer: 3,
    explanation: "Lactic and acetic acids and ethanol contribute to sourdough characteristics."
  },
  {
    question: "Biochemical oxygen demand (BOD) measures:",
    options: ["Oxygen consumed by microbes", "pH change", "COD removal", "Sludge volume"],
    answer: 0,
    explanation: "BOD estimates organic matter degraded by microbes consuming oxygen."
  },
  {
    question: "The enzyme cellulase is industrially valuable because it:",
    options: ["Hydrolyzes cellulose into sugars", "Synthesizes cellulose", "Denatures proteins", "Produces ethanol directly"],
    answer: 0,
    explanation: "Cellulases break down cellulose to fermentable sugars."
  },
  {
    question: "Which of the following microbes is a strict anaerobe important in digester methane production?",
    options: ["Methanobacterium", "Pseudomonas", "Bacillus", "Lactobacillus"],
    answer: 0,
    explanation: "Methanobacterium are methanogenic archaea thriving in strict anaerobic digesters."
  },
  {
    question: "Which method is used for rapid enumeration of viable bacteria in environmental samples?",
    options: ["qPCR", "Plate count", "Staining only", "Microscopy"],
    answer: 1,
    explanation: "Plate counts measure viable cells capable of colony formation."
  },
  {
    question: "Which process in biogeochemical cycling converts organic phosphorus to inorganic phosphate?",
    options: ["Mineralization", "Denitrification", "Nitrification", "Photosynthesis"],
    answer: 0,
    explanation: "Mineralization releases phosphate from organic compounds."
  },
  {
    question: "One major challenge of large‑scale algal biofuel production is:",
    options: ["Contamination control", "Maintaining sterile conditions", "Algal strain stability", "All of the above"],
    answer: 3,
    explanation: "Large‑scale cultivation faces contamination, strain drift, and stability issues."
  },
  {
    question: "Which bacterial group is best known for high‑yield amino acid production (e.g., glutamate)?",
    options: ["Corynebacterium", "Clostridium", "Enterococcus", "Vibrio"],
    answer: 0,
    explanation: "Corynebacterium glutamicum is industrially significant for amino acid fermentation."
  },
  {
    question: "In bioremediation of chlorinated solvents, which microbial process is essential?",
    options: ["Reductive dechlorination", "Oxidative phosphorylation", "Photosynthesis", "Denitrification"],
    answer: 0,
    explanation: "Reductive dechlorination removes chlorine atoms from solvents anaerobically."
  },
  {
    question: "Which microbial product commonly enhances soil nutrient availability?",
    options: ["Organic acids", "Aflatoxins", "Ammonia only", "H2S"],
    answer: 0,
    explanation: "Organic acids from microbes solubilize minerals, aiding plant uptake."
  },
  {
    question: "The Monod equation describes:",
    options: ["Microbial growth rate vs substrate", "Temperature effects", "pH control", "Foam formation"],
    answer: 0,
    explanation: "Monod kinetics relate growth rate to limiting substrate concentration."
  },
  {
    question: "Which method is used to sterilize large volumes in industry?",
    options: ["Autoclave only", "Continuous pasteurization", "Gamma irradiation", "Ozone sparging"],
    answer: 1,
    explanation: "Continuous pasteurization sterilizes flowing liquids at industrial scale."
  },
  {
    question: "In environmental microbiology, indicator microbes signal:",
    options: ["Specific pathogens", "Possible contamination", "Organic matter only", "Temperature change"],
    answer: 1,
    explanation: "Indicators like coliforms suggest fecal or pathogenic contamination."
  },
  {
    question: "Which of these is true about nitrification?",
    options: [
      "Occurs anaerobically",
      "Consumes oxygen",
      "Produces methane",
      "Produces alcohol"
    ],
    answer: 1,
    explanation: "Nitrifiers are aerobic and consume oxygen while oxidizing ammonia."
  },
  {
    question: "Which is NOT a function of microbial exoenzymes in industrial processes?",
    options: ["Protein degradation", "Cellulolysis", "DNA replication", "Starch hydrolysis"],
    answer: 2,
    explanation: "Exoenzymes act extracellularly, not for DNA replication."
  },
  {
    question: "Which method separates microbes based on density gradients?",
    options: ["Centrifugation", "PCR", "Gel electrophoresis", "Flow cytometry"],
    answer: 0,
    explanation: "Density gradient centrifugation separates cells by buoyant density."
  },
  {
    question: "Phosphate solubilizing microbes enhance plant growth by:",
    options: ["Fixing N2", "Releasing inorganic phosphate", "Producing H2S", "Reducing oxygen"],
    answer: 1,
    explanation: "They convert unavailable phosphate to forms plants can absorb."
  },
  {
    question: "Which organism is commonly used in industrial production of gluconic acid?",
    options: ["Aspergillus niger", "Lactobacillus acidophilus", "E. coli", "Saccharomyces cerevisiae"],
    answer: 0,
    explanation: "A. niger oxidizes glucose to gluconic acid at industrial scale."
  },
  {
    question: "Selective pressure in antibiotic fermentations is used to:",
    options: [
      "Increase yield", 
      "Prevent contamination", 
      "Maintain plasmid expression", 
      "All of the above"
    ],
    answer: 3,
    explanation: "Selective pressure helps maintain engineered traits and reduce contaminants."
  }
];
let quizData3 = [
  {
    question: "In the Embden–Meyerhof–Parnas (EMP) glycolytic pathway, which reaction is the primary irreversible, rate-limiting control point in most bacteria?",
    options: [
      "Hexokinase (glucose → glucose-6-phosphate)",
      "Phosphofructokinase-1 (fructose-6-phosphate → fructose-1,6-bisphosphate)",
      "Pyruvate kinase (phosphoenolpyruvate → pyruvate)",
      "Glyceraldehyde-3-phosphate dehydrogenase (G3P → 1,3-bisphosphoglycerate)"
    ],
    answer: 1,
    explanation: "Phosphofructokinase-1 (PFK-1) is the major committed, highly regulated step in EMP glycolysis. It is allosterically activated/inhibited by cellular energy signals (e.g., ADP/AMP activate; ATP and citrate inhibit), making it the central flux control point."
  },

  {
    question: "Which statement correctly explains why the Entner–Doudoroff (ED) pathway produces less net ATP per glucose than EMP?",
    options: [
      "ED bypasses substrate-level phosphorylation completely",
      "ED yields one NADPH and one NADH and generates only one net ATP per glucose",
      "ED oxidizes glucose with direct electron transfer to oxygen, producing less ATP",
      "ED produces three pyruvate molecules per glucose thereby lowering ATP yield"
    ],
    answer: 1,
    explanation: "The ED pathway yields one NADPH, one NADH and only one ATP (net) per glucose because it converts glucose to 6-phosphogluconate then to KDPG which cleaves to pyruvate and G3P; G3P proceeds to pyruvate with one ATP produced, so net ATP is lower than EMP's two ATP."
  },

  {
    question: "Which enzyme in the pentose phosphate pathway is the principal source of NADPH in bacteria for biosynthesis and oxidative stress defense?",
    options: [
      "Transketolase",
      "6-Phosphogluconate dehydrogenase",
      "Glucose-6-phosphate dehydrogenase (G6PD)",
      "Ribulose-5-phosphate epimerase"
    ],
    answer: 2,
    explanation: "Glucose-6-phosphate dehydrogenase catalyzes the first, rate-limiting oxidative step of the oxidative pentose phosphate pathway, producing NADPH that is used for biosynthesis (fatty acids, nucleotides) and for maintaining redox balance (e.g., glutathione reduction)."
  },

  {
    question: "In bacterial oxidative phosphorylation the majority of ATP synthesis occurs by which mechanism?",
    options: [
      "Substrate-level phosphorylation catalysed by kinases",
      "Direct oxidation of NADH to produce ATP",
      "Chemiosmotic coupling where the proton motive force drives ATP synthase",
      "Electron bifurcation with direct ATP-synthase activation"
    ],
    answer: 2,
    explanation: "Chemiosmotic coupling (Mitchell's chemiosmotic hypothesis) describes how electron transport pumps protons to create a proton motive force (PMF) across the membrane; ATP synthase uses PMF to synthesize ATP from ADP and Pi — the dominant mechanism in bacteria."
  },

  {
    question: "Which of the following terminal electron acceptors yields the least free-energy change (ΔG) when used in anaerobic respiration compared to O₂?",
    options: [
      "NO₃⁻ → N₂ (denitrification)",
      "SO₄²⁻ → H₂S (sulfate reduction)",
      "CO₂ → CH₄ (methanogenesis)",
      "Fe³⁺ → Fe²⁺ (iron reduction)"
    ],
    answer: 2,
    explanation: "Reduction of CO₂ to CH₄ (methanogenesis) is among the lowest-energy terminal electron accepting reactions listed, providing less free energy per electron transferred than denitrification or sulfate reduction. Methanogenesis is generally restricted to specialized archaea in strongly reducing environments."
  },

  {
    question: "Which of the following features of bacterial ATP synthase (F₁F₀ ATPase) is essential for coupling proton flow to ATP synthesis?",
    options: [
      "The c-ring rotation in F₀ and conformational changes transmitted to F₁",
      "Direct transfer of Pi through the F₀ channel",
      "ATP binding to membrane subunits in F₀",
      "Protonation of the β-subunit catalytic site"
    ],
    answer: 0,
    explanation: "Proton flow through F₀ causes rotation of the c-ring and central stalk (γ subunit) which induces conformational changes in the F₁ β subunits, coupling mechanical rotation to ATP synthesis. Protons do not directly enter the catalytic site; energy is transmitted mechanically."
  },

  {
    question: "Which regulator mediates the stringent response in bacteria by synthesizing (p)ppGpp under amino-acid starvation?",
    options: [
      "Sigma factor σ70",
      "RelA/SpoT homologs (RSH proteins)",
      "cAMP-CRP complex",
      "LexA repressor"
    ],
    answer: 1,
    explanation: "RelA and SpoT (and RelA/SpoT homologs) synthesize and hydrolyze the alarmone (p)ppGpp in response to nutrient stress (e.g., uncharged tRNA during amino-acid starvation). (p)ppGpp reprograms transcription, downregulating rRNA/tRNA synthesis and adjusting metabolism."
  },

  {
    question: "Which kinetic parameter indicates the substrate concentration at which the enzyme's reaction rate is half its maximum velocity and therefore relates to enzyme affinity?",
    options: [
      "Vmax",
      "Turnover number (kcat)",
      "Km",
      "kcat/Km"
    ],
    answer: 2,
    explanation: "Km (Michaelis constant) is the substrate concentration at which the reaction rate is half of Vmax and is often interpreted as an inverse measure of affinity (lower Km = higher affinity) under Michaelis-Menten assumptions."
  },

  {
    question: "In facultative anaerobes, which pair of enzymes is most important for detoxifying reactive oxygen species formed when cells encounter oxygen?",
    options: [
      "Catalase and superoxide dismutase (SOD)",
      "Periplasmic nitrate reductase and nitrite reductase",
      "Glutamine synthetase and nitrogenase",
      "Glyoxylate synthase and malate synthase"
    ],
    answer: 0,
    explanation: "Superoxide dismutase converts superoxide radical (O₂⁻) to H₂O₂, and catalase (or peroxidases) converts H₂O₂ to H₂O and O₂. These enzymes are critical for survival when facultative anaerobes are exposed to oxygen."
  },

  {
    question: "Which of the following best describes substrate-level phosphorylation during glycolysis?",
    options: [
      "ATP synthesis driven by proton gradient across membrane",
      "ATP generated when phosphate is transferred directly from a phosphorylated intermediate to ADP",
      "ATP produced by photophosphorylation",
      "ATP produced by the action of ATP synthase using sodium motive force"
    ],
    answer: 1,
    explanation: "Substrate-level phosphorylation refers to the direct transfer of a phosphate group from a phosphorylated metabolic intermediate (e.g., 1,3-bisphosphoglycerate or PEP) to ADP to form ATP, independent of membrane gradients."
  },

  {
    question: "Which transport system in bacteria couples the import of sugars to phosphotransfer that modifies the sugar at entry and integrates regulation of carbon metabolism?",
    options: [
      "ABC transporter",
      "Group translocation — Phosphoenolpyruvate-dependent phosphotransferase system (PTS)",
      "Symporter driven by Na⁺ gradient",
      "Type III secretion system"
    ],
    answer: 1,
    explanation: "The bacterial PTS uses phosphoenolpyruvate (PEP) as the phosphoryl donor to phosphorylate sugars during transport across the membrane. PTS both transports and chemically modifies the sugar and provides regulatory input to central metabolism (e.g., via EIIA phosphorylation state)."
  },

  {
    question: "Which compound is the immediate electron donor for nitrogenase during biological nitrogen fixation in diazotrophs?",
    options: [
      "NADH directly",
      "Ferredoxin (reduced)",
      "FADH₂",
      "ATP"
    ],
    answer: 1,
    explanation: "Reduced ferredoxin is the physiological electron donor to nitrogenase complex; electrons are funneled from ferredoxin (or flavodoxin) through reductase components to the nitrogenase FeMo-cofactor. ATP is required to drive electron transfer, but it is not the direct electron donor."
  },

  {
    question: "Which bacterial metabolic strategy best describes chemolithoautotrophy?",
    options: [
      "Use of organic carbon and organic electron donors for energy",
      "Use of light for energy and organic carbon for growth",
      "Oxidation of inorganic compounds for energy and fixation of CO₂ as carbon source",
      "Fermentation of sugars into organic acids"
    ],
    answer: 2,
    explanation: "Chemolithoautotrophs oxidize inorganic electron donors (e.g., H₂, NH₃, Fe²⁺) to generate energy and use CO₂ as their carbon source via autotrophic CO₂ fixation pathways (e.g., Calvin cycle, reverse TCA)."
  },

  {
    question: "Which mechanism most directly reduces the efficacy of many antibiotics within biofilms?",
    options: [
      "Random mutation rates are lower in biofilms",
      "EPS matrix and slow growth create diffusion barriers and phenotypic tolerance",
      "Biofilms exclusively consist of Gram-positive bacteria resistant to antibiotics",
      "Antibiotics are enzymatically inactivated only outside biofilms"
    ],
    answer: 1,
    explanation: "The extracellular polymeric substance (EPS) matrix of biofilms impedes antibiotic penetration, and cells in biofilms often have reduced growth/metabolic rates and can enter tolerant persister states — collectively causing decreased antibiotic efficacy (phenotypic tolerance) rather than only genetic resistance."
  },

  {
    question: "Which bacterial enzyme catalyzes the formation of peptidoglycan cross-links and is the primary target of β-lactam antibiotics?",
    options: [
      "Transglycosylase",
      "MurA",
      "Transpeptidase (penicillin-binding proteins, PBPs)",
      "D-alanine racemase"
    ],
    answer: 2,
    explanation: "Transpeptidases (a class of penicillin-binding proteins, PBPs) catalyze cross-linking of peptidoglycan peptide side chains. β-lactam antibiotics mimic the D-Ala-D-Ala substrate and irreversibly acylate PBPs, inhibiting transpeptidation and weakening the cell wall."
  },

  {
    question: "Which statement about bacterial two-component regulatory systems is most accurate?",
    options: [
      "Sensor kinase phosphorylates response regulator, which commonly binds DNA to alter transcription",
      "Response regulators are membrane-bound proteins that transport molecules",
      "Two-component systems always involve cAMP as second messenger",
      "They directly synthesize small RNAs to regulate translation"
    ],
    answer: 0,
    explanation: "Two-component systems consist of a membrane-associated sensor histidine kinase that detects environmental stimuli and autophosphorylates; it transfers phosphate to a cytoplasmic response regulator that often acts as a transcription factor to change gene expression."
  },

  {
    question: "Which metabolic adaptation allows some bacteria to grow at low water activity (high osmolarity) without plasmolysis?",
    options: [
      "Production of endospores",
      "Synthesis/accumulation of compatible solutes (osmoprotectants) such as proline, glycine betaine",
      "Increased membrane permeability to water",
      "Switch to anaerobic respiration only"
    ],
    answer: 1,
    explanation: "Bacteria facing high osmotic stress accumulate compatible solutes (small organic osmolytes) which balance internal osmotic pressure without interfering with cellular processes, preventing water loss and plasmolysis."
  },

  {
    question: "Which enzyme system in bacteria repairs UV-induced thymine dimers using a light-dependent mechanism?",
    options: [
      "Nucleotide excision repair (UvrABC)",
      "Photolyase (photo-reactivation)",
      "Mismatch repair (MutS/MutL)",
      "Base excision repair (glycosylases)"
    ],
    answer: 1,
    explanation: "Photolyase binds to UV-induced cyclobutane pyrimidine dimers and uses visible/near-UV light energy to cleave the dimer (photo-reactivation). Other systems such as NER can repair dimers in a light-independent manner."
  },

  {
    question: "Which of the following best explains the role of alternative sigma factors in bacterial transcription?",
    options: [
      "They enhance ribosome assembly for increased translation",
      "They bind specific promoter sequences to enable expression of genes for stress or developmental programs",
      "They are proteases that degrade misfolded proteins during stress",
      "They serve as chaperones for DNA replication"
    ],
    answer: 1,
    explanation: "Alternative sigma factors replace the primary sigma factor in the RNA polymerase holoenzyme to redirect transcription to specific sets of promoters, allowing rapid reprogramming of gene expression during stress (e.g., σS for stationary phase, σH for heat shock)."
  },

  {
    question: "Which biochemical change is most characteristic of the bacterial stationary phase compared to exponential growth?",
    options: [
      "Increased ribosomal RNA transcription",
      "Upregulation of stress-response proteins and accumulation of (p)ppGpp",
      "Maximal rates of DNA replication and protein synthesis",
      "Complete shutdown of glycolysis"
    ],
    answer: 1,
    explanation: "During stationary phase, nutrient limitation triggers stress responses including increased (p)ppGpp (stringent response), downregulation of rRNA synthesis, and upregulation of stress-protective proteins (chaperones, antioxidant enzymes), enabling survival in non-growing states."
  },

  {
    question: "Which metabolic route is primarily used by many bacteria for anabolic (biosynthetic) generation of sugars and precursors from CO₂?",
    options: [
      "Glyoxylate shunt",
      "Calvin-Benson-Bassham cycle (ribulose bisphosphate carboxylase/oxygenase (RuBisCO)-dependent fixation)",
      "Embden–Meyerhof–Parnas pathway in reverse",
      "Entner–Doudoroff pathway running backward"
    ],
    answer: 1,
    explanation: "The Calvin cycle, employing RuBisCO and phosphoribulokinase, is the major CO₂-fixation pathway in many autotrophic bacteria and cyanobacteria, producing triose phosphates used for anabolic processes; other autotrophic pathways exist but Calvin is primary."
  },

  {
    question: "Why do anaerobic respiration processes (e.g., using nitrate or sulfate) generally yield less ATP than aerobic respiration using O₂?",
    options: [
      "Anaerobic respiration uses substrate-level phosphorylation exclusively",
      "Non-oxygen terminal electron acceptors have lower redox potential differences with donors, providing less free energy for proton pumping",
      "Electron transport chains are absent in anaerobes",
      "Anaerobic respiring organisms lack ATP synthase"
    ],
    answer: 1,
    explanation: "The free-energy change available from electron transfer depends on the redox potential difference between donor and acceptor. O₂ has a very high positive redox potential; many alternative acceptors (NO₃⁻, SO₄²⁻) have lower potentials, yielding less energy per electron pair and therefore fewer protons pumped and less ATP."
  },

  {
    question: "Which factor most directly determines the direction of proton-driven flagellar rotation and thus chemotactic behavior in bacteria like E. coli?",
    options: [
      "ATP concentration in the cytoplasm",
      "The phosphorylation state of chemotaxis proteins (e.g., CheY-P interacting with the flagellar motor)",
      "Presence of pili",
      "Cell membrane thickness"
    ],
    answer: 1,
    explanation: "CheY when phosphorylated (CheY-P) binds to the flagellar motor switch complex and changes rotation bias (switching between run and tumble). Proton flow powers rotation, but switching direction for chemotaxis is controlled by chemotaxis signaling proteins like CheY."
  },

  {
    question: "In the context of catabolite repression in E. coli, what role does cAMP-CRP play for transcription of alternative carbon-source operons (e.g., lac operon)?",
    options: [
      "cAMP-CRP complex inhibits RNA polymerase binding to promoters",
      "Low glucose increases cAMP; cAMP-CRP binds promoters to activate transcription of catabolic operons",
      "High glucose increases cAMP levels and activates catabolic genes",
      "cAMP-CRP cleaves operator DNA to permit transcription"
    ],
    answer: 1,
    explanation: "When glucose is scarce, intracellular cAMP rises; cAMP forms a complex with CRP (CAP) which binds promoter-proximal sites and facilitates RNA polymerase recruitment, enhancing transcription of operons needed to utilize alternative carbon sources. High glucose reduces cAMP, repressing these operons (catabolite repression)."
  },

  {
    question: "Which enzyme is directly responsible for initiating peptidoglycan monomer synthesis in the cytoplasm by catalyzing the first committed step of the UDP-N-acetylglucosamine pathway?",
    options: [
      "MurA (enolpyruvyl transferase)",
      "Transglycosylase",
      "GlmS (glutamine—fructose-6-phosphate aminotransferase)",
      "D-Ala-D-Ala ligase (Ddl)"
    ],
    answer: 0,
    explanation: "MurA (UDP-N-acetylglucosamine enolpyruvyl transferase) catalyzes the transfer of enolpyruvate from phosphoenolpyruvate to UDP-N-acetylglucosamine, committing the substrate toward peptidoglycan synthesis. It is targeted by the antibiotic fosfomycin."
  }

];

let quizData4 = [
  {
    question: "Which protein in E. coli recognizes the origin of replication (oriC) and initiates assembly of the replication fork?",
    options: [
      "DnaA",
      "DnaB (helicase)",
      "DNA polymerase III",
      "Primase (DnaG)"
    ],
    answer: 0,
    explanation: "DnaA is the origin-binding initiator protein in E. coli. It binds DnaA-boxes at oriC, oligomerizes and melts the AT-rich region to recruit the helicase (DnaB) and other replisome components, initiating bidirectional replication."
  },
  {
    question: "During bacterial DNA replication, which enzyme is primarily responsible for elongation on the leading and lagging strands?",
    options: [
      "DNA polymerase I",
      "DNA polymerase II",
      "DNA polymerase III holoenzyme",
      "Primase"
    ],
    answer: 2,
    explanation: "DNA polymerase III holoenzyme is the main replicative polymerase in bacteria; it synthesizes both the continuous leading strand and the Okazaki fragments on the lagging strand with high processivity (clamp) and proofreading (ε subunit). DNA pol I removes RNA primers and fills gaps."
  },
  {
    question: "Which molecular event directly generates the discontinuous nature of the lagging strand during replication?",
    options: [
      "RNA primer synthesis by primase",
      "Helicase translocation",
      "DNA gyrase relieving supercoils",
      "DNA ligase sealing nicks"
    ],
    answer: 0,
    explanation: "Primase synthesizes short RNA primers at intervals on the lagging-strand template; DNA pol III extends from these primers to form Okazaki fragments, producing discontinuous synthesis. Subsequent removal of primers and ligation yields a continuous strand."
  },
  {
    question: "In conjugative transfer mediated by an F plasmid, what is the defining feature of an Hfr strain?",
    options: [
      "It lacks the tra genes",
      "It carries the F factor integrated into the chromosome",
      "It is naturally competent for transformation",
      "It produces generalized transducing phage"
    ],
    answer: 1,
    explanation: "Hfr (high frequency recombination) strains have the F plasmid integrated into the host chromosome. During conjugation they can transfer chromosomal genes beginning at the integrated oriT, enabling transfer of donor chromosomal loci to recipients."
  },
  {
    question: "Which type of transduction transfers only genes adjacent to the prophage integration site and is dependent on improper excision of the prophage?",
    options: [
      "Generalized transduction",
      "Specialized transduction",
      "Conjugative transduction",
      "Transformational transduction"
    ],
    answer: 1,
    explanation: "Specialized transduction occurs when a temperate phage (prophage) excises imprecisely and carries specific host genes adjacent to the prophage insertion site into new host cells. Generalized transduction packages random host DNA during lytic cycles."
  },
  {
    question: "Which bacterial system mediates homologous recombination and is essential for strand invasion during repair and recombination?",
    options: [
      "RecA",
      "MutS",
      "UvrABC",
      "LigA"
    ],
    answer: 0,
    explanation: "RecA polymerizes on single-stranded DNA and catalyzes strand invasion and exchange during homologous recombination; it also facilitates autocleavage of the LexA repressor to induce the SOS response under DNA damage."
  },
  {
    question: "Which enzyme complex removes RNA primers and replaces them with DNA during bacterial replication?",
    options: [
      "DNA polymerase III",
      "DNA polymerase I",
      "RNase H alone",
      "DNA ligase"
    ],
    answer: 1,
    explanation: "DNA polymerase I has 5'→3' exonuclease activity that removes RNA primers and 5'→3' polymerase activity that fills the resulting gaps with DNA. DNA ligase then seals the remaining nicks."
  },
  {
    question: "Riboswitches regulate gene expression by binding small metabolites. At which level do most riboswitches act in bacteria?",
    options: [
      "DNA replication initiation",
      "Transcription termination or translation initiation",
      "Protein folding",
      "Protein degradation"
    ],
    answer: 1,
    explanation: "Riboswitches are RNA elements usually located in the 5' UTR of mRNAs; binding of a metabolite causes a structural change that affects transcription termination (attenuation) or blocks/unblocks the ribosome binding site, thereby regulating gene expression at the transcriptional or translational level."
  },
  {
    question: "Which enzyme is the direct target of the antibiotic rifampicin, and why does it inhibit bacterial transcription?",
    options: [
      "Sigma factor; prevents promoter recognition",
      "RNA polymerase β subunit; blocks transcript elongation",
      "DNA gyrase; prevents supercoiling",
      "RNA helicase; prevents initiation complex formation"
    ],
    answer: 1,
    explanation: "Rifampicin binds to the β subunit of bacterial RNA polymerase and sterically blocks the path of the nascent RNA chain, preventing transcription elongation after initiation. It does not block promoter binding but inhibits early RNA synthesis."
  },
  {
    question: "Which mechanism best describes attenuation control as exemplified by the trp operon in E. coli?",
    options: [
      "A protein repressor binding operator to stop RNA polymerase",
      "Ribosome-mediated modulation of transcription termination in response to charged tRNA levels",
      "DNA methylation preventing transcription factor binding",
      "Small RNA binding to mRNA to prevent translation"
    ],
    answer: 1,
    explanation: "Attenuation in the trp operon uses translation of a leader peptide; when tryptophan is abundant, ribosomes translate the leader quickly, enabling formation of a rho-independent terminator hairpin in the mRNA that halts transcription. When tryptophan is scarce, ribosome stalling alters RNA structure to allow transcription of downstream genes."
  },
  {
    question: "CRISPR-Cas adaptive immunity in bacteria stores short sequences of invading DNA in the CRISPR array. Which step uses the stored spacers to degrade incoming matching nucleic acids?",
    options: [
      "Adaptation",
      "Expression and processing into crRNAs followed by interference",
      "Insertion of spacers into host genome",
      "Conjugation"
    ],
    answer: 1,
    explanation: "After adaptation (spacer acquisition), CRISPR arrays are transcribed and processed into CRISPR RNAs (crRNAs). During interference, crRNAs guide Cas nucleases to complementary foreign nucleic acid sequences for cleavage, providing sequence-specific immunity."
  },
  {
    question: "Which polymerase is typically used in standard PCR because of its thermostability and lack of 3'→5' exonuclease activity sufficient for high-temperature cycling?",
    options: [
      "E. coli DNA polymerase I",
      "Taq DNA polymerase",
      "DNA polymerase III",
      "Reverse transcriptase"
    ],
    answer: 1,
    explanation: "Taq DNA polymerase, isolated from Thermus aquaticus, is thermostable and functions at high temperatures used in PCR. Many PCR variants use proofreading polymerases (with 3'→5' exonuclease activity) for higher fidelity, but standard PCR historically used Taq."
  },
  {
    question: "Which feature most characterizes rolling-circle replication used by many plasmids and some phages?",
    options: [
      "Bidirectional forks starting at oriC",
      "A single-strand nick that provides a 3'-OH primer and continuous synthesis displacing the old strand",
      "Replication entirely by host primase without helicase",
      "Synthesis of Okazaki fragments on both strands"
    ],
    answer: 1,
    explanation: "Rolling-circle replication begins with a site-specific nick in one strand; DNA synthesis extends the 3' end displacing the old strand, producing a single-stranded linear DNA that can be converted to double-stranded form. It is unidirectional and continuous for the leading strand."
  },
  {
    question: "Which mutational change is most likely to have the largest structural impact on a protein and frequently results from insertion or deletion of nucleotides not divisible by three?",
    options: [
      "Missense mutation",
      "Nonsense mutation",
      "Frameshift mutation",
      "Silent mutation"
    ],
    answer: 2,
    explanation: "Frameshift mutations caused by insertions or deletions that are not multiples of three alter the reading frame downstream, changing all subsequent codons and usually producing a nonfunctional protein or early stop codon with dramatic structural consequences."
  },
  {
    question: "Which bacterial DNA repair pathway directly recognizes and removes bulky helix-distorting lesions such as thymine dimers in a light-independent manner?",
    options: [
      "Photoreactivation",
      "Nucleotide excision repair (UvrABC excinuclease)",
      "Mismatch repair (MutS/MutL)",
      "Base excision repair (glycosylases)"
    ],
    answer: 1,
    explanation: "Nucleotide excision repair (NER), involving UvrA/UvrB/UvrC in bacteria, removes bulky, helix-distorting lesions (including UV-induced thymine dimers) by excising an oligonucleotide segment around the lesion; DNA polymerase and ligase fill and seal the gap."
  },
  {
    question: "In plasmid biology, which region contains the minimal elements required for autonomous replication in a given host?",
    options: [
      "tra operon",
      "oriT",
      "origin of replication (ori) and replication initiator gene(s)",
      "Antibiotic resistance gene"
    ],
    answer: 2,
    explanation: "Autonomous plasmid replication requires an origin of replication (ori) and often plasmid-encoded initiator proteins (Rep proteins) that recognize ori and recruit host replication machinery. tra genes mediate transfer but are not required for replication per se."
  },
  {
    question: "What is the principal molecular basis for codon bias observed between highly expressed and lowly expressed genes in bacteria?",
    options: [
      "Random mutation of codons",
      "Selection for codons that match abundant tRNAs to increase translation efficiency and accuracy",
      "tRNA degradation mechanisms",
      "Ribosomal RNA processing differences"
    ],
    answer: 1,
    explanation: "Codon bias arises because highly expressed genes preferentially use synonymous codons that correspond to abundant tRNA isoacceptors in the cell, enhancing translation speed and fidelity; selection for translational efficiency shapes codon usage patterns."
  },
  {
    question: "Which small, non-coding RNAs (sRNAs) commonly regulate bacterial mRNA stability or translation by base-pairing and often require the Hfq protein as a chaperone?",
    options: [
      "Ribosomal RNA fragments",
      "CRISPR RNAs",
      "Trans-acting sRNAs (e.g., RyhB) that modulate target mRNAs",
      "mRNA leader peptides"
    ],
    answer: 2,
    explanation: "Trans-acting sRNAs are small regulatory RNAs that base-pair with target mRNAs to influence translation initiation or stability; Hfq often stabilizes sRNA-mRNA interactions. Examples include RyhB (iron metabolism) in E. coli."
  },
  {
    question: "Which viral genome replication strategy creates double-stranded DNA intermediates from single-stranded DNA viruses for transcription by host RNA polymerase?",
    options: [
      "Positive-sense ssRNA replication",
      "Rolling-circle replication of ssDNA viruses producing dsDNA replicative forms",
      "Retroviral reverse transcription into dsDNA",
      "Negative-sense ssRNA replication without DNA intermediates"
    ],
    answer: 1,
    explanation: "Many ssDNA viruses (e.g., parvoviruses, some bacteriophages) use rolling-circle or strand displacement mechanisms to create double-stranded replicative forms that serve as templates for transcription by host RNA polymerase; these dsDNA intermediates are necessary for mRNA synthesis."
  },
  {
    question: "Which molecular technique allows precise editing of bacterial chromosomes by creating a double-strand break at a chosen site followed by homologous recombination with an engineered donor template?",
    options: [
      "Random transposon mutagenesis",
      "CRISPR-Cas9–directed cleavage combined with homologous recombination (recombineering)",
      "Blue-white screening",
      "Southern blotting"
    ],
    answer: 1,
    explanation: "CRISPR-Cas9 can be programmed to introduce a targeted double-strand break; providing a donor DNA with homology arms enables homologous recombination-based repair (recombineering), allowing precise genome edits. This strategy is widely used for bacterial genome engineering."
  },
  {
    question: "Which conserved sequence element in bacterial promoters is typically recognized by the σ70 family sigma factor at about 10 bases upstream of the transcription start site?",
    options: [
      "The −35 element (TTGACA)",
      "The Shine-Dalgarno sequence",
      "The −10 (Pribnow) box (TATAAT)",
      "Terminator hairpin"
    ],
    answer: 2,
    explanation: "The −10 element (Pribnow box, consensus TATAAT) is recognized along with the −35 element by σ70-containing RNA polymerase holoenzyme and is typically located ~10 bases upstream of the transcription start site; together these elements constitute core promoter recognition motifs."
  },
  {
    question: "Which genomic feature commonly carried on integrons contributes to rapid spread of antibiotic resistance among bacteria?",
    options: [
      "Conserved ribosomal RNA operons",
      "Gene cassettes and site-specific integrase that capture and express resistance genes",
      "Replication origins",
      "CRISPR arrays"
    ],
    answer: 1,
    explanation: "Integrons are genetic platforms that capture mobile gene cassettes via a site-specific integrase (IntI); many cassettes carry antibiotic resistance genes which, when expressed from the integron promoter, facilitate rapid acquisition and dissemination of resistance."
  },
  {
    question: "What is the primary role of the LexA protein in bacterial cells under normal (non-DNA-damage) conditions?",
    options: [
      "Activates SOS genes",
      "Represses SOS response genes by binding operator sequences",
      "Functions as a DNA helicase during replication",
      "Binds single-stranded DNA to initiate recombination"
    ],
    answer: 1,
    explanation: "LexA is a repressor that binds SOS box sequences in promoters of DNA-damage-inducible genes, repressing their transcription. Upon DNA damage, RecA stimulates autocleavage of LexA, derepressing SOS genes involved in repair and tolerance."
  },
  {
    question: "Which sequencing approach is best suited for detecting structural genome variations (large insertions, deletions, inversions) and resolving repetitive regions in microbial genomes?",
    options: [
      "Short-read sequencing (e.g., Illumina) alone",
      "Long-read sequencing technologies (e.g., PacBio, Oxford Nanopore)",
      "Sanger sequencing of PCR amplicons only",
      "Microarray hybridization"
    ],
    answer: 1,
    explanation: "Long-read sequencing (PacBio, Oxford Nanopore) produces reads that span large repetitive elements and structural variants, facilitating assembly of complete genomes and detection of insertions, deletions, inversions and complex rearrangements that are difficult to resolve with short reads alone."
  }
];

let quizData5 = [
  {
    question: "Which intracellular signaling event is the immediate result of TCR engagement with peptide-MHC and co-stimulation that leads to full T cell activation?",
    options: [
      "Activation of JAK-STAT pathway leading to gene transcription",
      "Phosphorylation of ITAMs on CD3 chains by Lck and subsequent ZAP-70 recruitment",
      "Direct increase of intracellular cAMP to activate PKA",
      "Cleavage of pro-IL-1β to IL-1β by caspase-1"
    ],
    answer: 1,
    explanation: "TCR engagement triggers Lck (a Src family kinase) to phosphorylate ITAM motifs on CD3 chains, which recruits and activates ZAP-70. ZAP-70 phosphorylates downstream adaptor proteins (LAT, SLP-76), initiating Ca2+, MAPK and NF-κB pathways required for full T cell activation; JAK-STAT is more typical of cytokine receptor signaling."
  },
  {
    question: "Which molecular event is required for Class I MHC presentation of viral peptides derived from cytosolic proteins?",
    options: [
      "Proteolytic cleavage in endosomes by cathepsins",
      "Proteasomal degradation, TAP transport into ER, and loading onto MHC I",
      "Processing in the Golgi with mannose trimming",
      "Autophagic delivery of peptides directly to MHC II"
    ],
    answer: 1,
    explanation: "Endogenous proteins (e.g., viral) are ubiquitinated and degraded by the proteasome; resulting peptides are transported into the ER by TAP, trimmed and loaded onto MHC class I molecules for presentation to CD8+ T cells. Endosomal cathepsins are used for MHC II."
  },
  {
    question: "Which cytokine is most critical for Th17 differentiation from naive CD4+ T cells in humans (in combination with TGF-β)?",
    options: [
      "IL-12",
      "IL-4",
      "IL-6 (and IL-23 for maintenance)",
      "IL-2"
    ],
    answer: 2,
    explanation: "In humans, naive CD4+ T cells differentiate into Th17 cells in the presence of TGF-β plus pro-inflammatory cytokines like IL-6; IL-23 is important for stabilization/maintenance of Th17 cells. IL-12 drives Th1, IL-4 drives Th2, IL-2 supports T cell proliferation."
  },
  {
    question: "Which complement regulator on host cells prevents formation of C3 convertase and protects self cells from complement-mediated damage?",
    options: [
      "Properdin (factor P)",
      "Factor B",
      "Decay-accelerating factor (DAF/CD55)",
      "C5a receptor"
    ],
    answer: 2,
    explanation: "DAF (CD55) expressed on host cell surfaces accelerates dissociation of C3 convertases (both classical/lectin and alternative pathway C3 convertases), preventing complement activation on self cells. Properdin stabilizes convertase, factor B is part of alternative pathway, and C5aR is a receptor for the anaphylatoxin C5a."
  },
  {
    question: "Affinity maturation in germinal centers primarily depends on which enzyme that introduces mutations into immunoglobulin variable regions?",
    options: [
      "Activation-induced cytidine deaminase (AID)",
      "RAG1/RAG2 recombinase",
      "Terminal deoxynucleotidyl transferase (TdT)",
      "DNA polymerase I"
    ],
    answer: 0,
    explanation: "AID deaminates cytosine to uracil in Ig variable region DNA during somatic hypermutation and also initiates class-switch recombination; RAG mediates V(D)J recombination in developing lymphocytes, TdT adds N nucleotides during V(D)J recombination."
  },
  {
    question: "Which innate immune receptor family senses cytosolic DNA and induces type I interferon production via STING?",
    options: [
      "Toll-like receptors (TLR1/2/6)",
      "NOD-like receptors (NOD1/NOD2)",
      "cGAS (cyclic GMP-AMP synthase) — STING pathway",
      "RIG-I-like receptors (RIG-I/MDA5)"
    ],
    answer: 2,
    explanation: "cGAS recognizes cytosolic double-stranded DNA and synthesizes cyclic GMP-AMP (cGAMP), which activates STING on the ER to trigger TBK1-IRF3 signaling and type I interferon production. RIG-I/MDA5 detect cytosolic RNA, TLRs detect extracellular/endosomal PAMPs, and NODs sense bacterial peptidoglycan fragments."
  },
  {
    question: "Which immunoglobulin is most effective at opsonization and complement activation in human serum?",
    options: [
      "IgA",
      "IgD",
      "IgG (particularly IgG1 and IgG3 subclasses)",
      "IgE"
    ],
    answer: 2,
    explanation: "IgG, especially subclasses IgG1 and IgG3, are efficient at opsonizing pathogens for phagocytosis and activating the classical complement pathway via C1q binding. IgA is important at mucosal surfaces, IgE in allergy/parasite defense, IgD has limited understood effector roles."
  },
  {
    question: "Which cell type is the principal professional antigen-presenting cell capable of priming naive CD8+ T cells via cross-presentation?",
    options: [
      "Neutrophils",
      "B cells",
      "Conventional dendritic cells (cDC1 subset)",
      "Eosinophils"
    ],
    answer: 2,
    explanation: "Certain dendritic cell subsets (notably cDC1 in mice/humans) are specialized for cross-presentation — presenting exogenous antigens on MHC I to prime naive CD8+ T cells; B cells and macrophages present antigen but are less efficient at cross-priming naive CD8+ T cells."
  },
  {
    question: "Which immune checkpoint receptor expressed on T cells transmits an inhibitory signal upon binding its ligand and is targeted by cancer immunotherapy (checkpoint blockade)?",
    options: [
      "CD28",
      "PD-1 (programmed cell death protein 1)",
      "CD40L",
      "IL-2R"
    ],
    answer: 1,
    explanation: "PD-1 is an inhibitory receptor on activated T cells; engagement by PD-L1/PD-L2 on tumor or other cells reduces T cell effector functions. Antibodies blocking PD-1/PD-L1 (checkpoint inhibitors) can restore antitumor T cell responses; CD28 is co-stimulatory."
  },
  {
    question: "In a primary immune response to a protein antigen, which antibody is produced first and predominates initially in serum?",
    options: [
      "IgG (high affinity, class-switched)",
      "IgA (mucosal dimeric form)",
      "IgM (pentameric, produced by naive B cells)",
      "IgE (associated with allergy)"
    ],
    answer: 2,
    explanation: "During a primary T-dependent response, naive B cells initially secrete IgM (low affinity, pentameric), which can efficiently activate complement. Over time and with T cell help, class switching produces IgG, IgA or IgE with higher affinity (after somatic hypermutation and selection)."
  },
  {
    question: "Which laboratory assay quantitatively measures cytokine secretion from single cells and allows enumeration of antigen-specific T cells by spot formation?",
    options: [
      "ELISA",
      "ELISPOT",
      "Western blot",
      "Flow cytometry scatter plots"
    ],
    answer: 1,
    explanation: "ELISPOT detects and counts individual cells secreting a specific cytokine (e.g., IFN-γ) by capturing secreted cytokine on a membrane and visualizing spots; it's highly sensitive for enumerating antigen-specific T cells. ELISA measures bulk cytokine concentration in supernatant."
  },
  {
    question: "Which mechanism explains central tolerance of developing T cells in the thymus to many tissue-specific antigens?",
    options: [
      "Expression of tissue-restricted antigens by medullary thymic epithelial cells (AIRE-dependent)",
      "Sequestration of antigens in peripheral immune privileged sites only",
      "Peripheral anergy induction by lack of co-stimulation",
      "Somatic hypermutation deletes autoreactive TCRs"
    ],
    answer: 0,
    explanation: "AIRE (autoimmune regulator) in medullary thymic epithelial cells drives promiscuous expression of many tissue-specific antigens, enabling negative selection (deletion) of T cells with high affinity for self. Peripheral mechanisms handle remaining autoreactive cells."
  },
  {
    question: "Which clinical hypersensitivity reaction is primarily mediated by immune complexes that deposit in tissues and activate complement leading to inflammation?",
    options: [
      "Type I hypersensitivity (IgE mediated)",
      "Type II hypersensitivity (antibody mediated cytotoxic)",
      "Type III hypersensitivity (immune complex mediated)",
      "Type IV hypersensitivity (cell mediated delayed-type)"
    ],
    answer: 2,
    explanation: "Type III hypersensitivity involves formation of antigen-antibody complexes that can deposit in vessel walls or tissues, activate complement and recruit neutrophils causing inflammation, as seen in serum sickness or some forms of glomerulonephritis."
  },
  {
    question: "Which signaling pathway is critical for B cell isotype switching (class switch recombination) and depends on activation-induced cytidine deaminase (AID)?",
    options: [
      "Somatic hypermutation via RAG recombinase",
      "Class switch recombination guided by cytokine milieu and AID-mediated DNA deamination",
      "Ig light-chain gene rearrangement by TdT",
      "V(D)J recombination in peripheral B cells"
    ],
    answer: 1,
    explanation: "CSR is initiated by AID deaminating cytosines in switch regions upstream of constant region genes; combined with double-strand break repair pathways and guidance from cytokines (e.g., IL-4 promotes IgG1/IgE switching), this results in a new heavy-chain constant region and altered effector function."
  },
  {
    question: "Which innate cell type kills virus-infected cells by detecting absence of MHC class I (missing-self) or stress ligands and does not require prior sensitization?",
    options: [
      "CD8+ cytotoxic T lymphocyte",
      "Natural killer (NK) cell",
      "Neutrophil",
      "B cell"
    ],
    answer: 1,
    explanation: "NK cells recognize cells with reduced MHC I expression (a common viral evasion tactic) or stress-induced ligands via activating receptors; they mediate cytotoxicity without prior antigen-specific sensitization and secrete IFN-γ to shape adaptive responses."
  },
  {
    question: "Which technique directly measures the diversity of a patient's T cell receptor (TCR) repertoire by sequencing and can reveal clonal expansions?",
    options: [
      "ELISA for cytokines",
      "TCR deep sequencing (immune repertoire sequencing)",
      "Immunohistochemistry for CD3",
      "Western blot for TCR beta"
    ],
    answer: 1,
    explanation: "High-throughput (deep) sequencing of TCR clonotypes (VDJ regions) quantitatively assesses T cell diversity and clonality, detecting expanded clones (e.g., in infection, cancer or autoimmune disease). ELISA measures soluble proteins, not receptor diversity."
  },
  {
    question: "Which adjuvant mechanism most commonly enhances vaccine immunogenicity by activating innate immune receptors and promoting stronger adaptive responses?",
    options: [
      "Direct activation of B cell antibody secretion without T cell help",
      "Depot effect only (slow antigen release) without innate activation",
      "Stimulation of pattern recognition receptors (e.g., TLR agonists) on APCs to upregulate co-stimulation and cytokines",
      "Blocking MHC presentation to prolong antigen presence"
    ],
    answer: 2,
    explanation: "Many adjuvants (e.g., TLR agonists like CpG) activate innate PRRs on antigen-presenting cells, increasing expression of co-stimulatory molecules and cytokines that enhance T and B cell priming; depot effects may contribute but PRR stimulation is key for robust adaptive responses."
  },
  {
    question: "Which of the following best describes antibody-dependent cellular cytotoxicity (ADCC)?",
    options: [
      "Complement-mediated lysis after antibody binding",
      "NK cells recognize Fc portions of antibodies bound to target via FcγRIII (CD16) and induce target cell apoptosis",
      "Direct neutralization of toxin by Fab binding",
      "Phagocytosis by macrophages via complement receptors only"
    ],
    answer: 1,
    explanation: "ADCC occurs when Fc receptor-bearing effector cells (e.g., NK cells via FcγRIII/CD16) bind the Fc region of antibodies coating target cells; this engagement triggers release of cytotoxic granules (perforin/granzymes) leading to apoptosis of the target. Complement-mediated lysis is distinct."
  },
  {
    question: "Which laboratory method is used to functionally assess neutrophil oxidative burst capacity and diagnose chronic granulomatous disease (CGD)?",
    options: [
      "Nitroblue tetrazolium (NBT) test or dihydrorhodamine (DHR) flow cytometry assay",
      "ELISA for neutrophil elastase",
      "Western blot for NADPH oxidase subunits",
      "Serum complement C3 measurement"
    ],
    answer: 0,
    explanation: "The NBT reduction test (historical) and DHR flow cytometry assay (modern) measure neutrophil respiratory burst by detecting reactive oxygen species production; defective oxidative burst indicates CGD due to NADPH oxidase defects. Western blot may detect subunits but functional assays are diagnostic."
  },
  {
    question: "Which mechanism allows peripheral tolerance mediated by regulatory T cells (Tregs) to suppress autoreactive immune responses?",
    options: [
      "Secretion of anti-inflammatory cytokines (IL-10, TGF-β), CTLA-4 mediated modulation of APCs, and metabolic disruption",
      "Induction of hypermutation to eliminate autoreactive receptors",
      "Promotion of complement activation on self cells",
      "Enhancing DC maturation to boost immune responses"
    ],
    answer: 0,
    explanation: "Tregs suppress immunity via multiple mechanisms including secretion of IL-10/TGF-β, expression of CTLA-4 to downregulate costimulatory molecules on APCs, consumption of IL-2 (cytokine sink), and metabolic interference — all contributing to peripheral tolerance and prevention of autoimmunity."
  },
  {
    question: "Which serum finding is most indicative of a humoral primary immunodeficiency characterized by failure of B cell development (e.g., X-linked agammaglobulinemia)?",
    options: [
      "Normal levels of all immunoglobulin isotypes",
      "Markedly reduced/absent serum IgG, IgM, and IgA with very low B cell counts",
      "Elevated IgE only",
      "Isolated low IgA with normal IgG/IgM"
    ],
    answer: 1,
    explanation: "X-linked agammaglobulinemia (BTK deficiency) leads to failure of B cell maturation — patients have very low/absent peripheral B cells and markedly reduced all immunoglobulin isotypes (IgG, IgM, IgA), causing recurrent bacterial infections. Isolated IgA deficiency or elevated IgE are different conditions."
  },
  {
    question: "Which pathway of complement activation is initiated by antibody-antigen complexes and requires C1q binding to clustered IgG or IgM?",
    options: [
      "Alternative pathway",
      "Lectin pathway",
      "Classical pathway",
      "Properdin pathway"
    ],
    answer: 2,
    explanation: "The classical complement pathway is initiated when C1q binds to Fc regions of antigen-bound IgM or clustered IgG, leading to activation of C1r/C1s, cleavage of C4/C2 to form C3 convertase. The lectin pathway is initiated by MBL binding sugars, and the alternative pathway is spontaneously activated on non-self surfaces."
  },
  {
    question: "Which cell surface marker combination is typically used to identify human regulatory T cells (Tregs) by flow cytometry?",
    options: [
      "CD3+ CD8+ CD56+",
      "CD4+ CD25high FOXP3+ (intracellular)",
      "CD19+ CD27+",
      "CD4+ IL-17+"
    ],
    answer: 1,
    explanation: "Human Tregs are commonly identified as CD4+ T cells with high CD25 (IL-2Rα) expression and intracellular FOXP3 transcription factor; FOXP3 is the lineage-defining regulator for Tregs. CD19 marks B cells, CD8 marks cytotoxic T cells, IL-17 marks Th17 cells."
  },
  {
    question: "Which phenomenon best explains why some vaccines (e.g., conjugate polysaccharide vaccines) are more immunogenic in infants than plain polysaccharide vaccines?",
    options: [
      "Conjugation converts a T-independent antigen into a T-dependent one by linking polysaccharide to protein carrier allowing T cell help and memory",
      "Conjugation masks the antigen preventing immune recognition",
      "Conjugated vaccines stimulate innate immune receptors only",
      "Polysaccharide vaccines are always superior in infants"
    ],
    answer: 0,
    explanation: "Polysaccharide antigens alone often elicit T-independent B cell responses with poor memory and weak responses in infants. Conjugate vaccines link the polysaccharide to a protein carrier, enabling processing and presentation to helper T cells and generation of class switching, affinity maturation and memory — effective in infants."
  }
];


let quizData6 = [
  {
    question: "Which virulence factor allows Staphylococcus aureus to evade opsonization and phagocytosis?",
    options: [
      "Protein A binding Fc portion of IgG",
      "Lipopolysaccharide endotoxin",
      "M protein",
      "Catalase enzyme"
    ],
    answer: 0,
    explanation: "Protein A binds the Fc region of IgG antibodies, preventing opsonization and phagocytosis by neutrophils and macrophages. LPS is a Gram-negative endotoxin, M protein is a virulence factor of Streptococcus pyogenes, and catalase breaks down hydrogen peroxide but is not directly involved in immune evasion."
  },
  {
    question: "Which secretion system in Gram-negative bacteria injects effector proteins directly into host cells to manipulate cellular processes?",
    options: [
      "Type I secretion system",
      "Type III secretion system",
      "Type V secretion system",
      "Sec pathway"
    ],
    answer: 1,
    explanation: "Type III secretion system acts as a molecular syringe to inject bacterial effector proteins into host cells, altering cytoskeleton, signaling, and immune responses. Type I and V systems secrete proteins extracellularly, and Sec pathway translocates proteins across the inner membrane but not directly into host cells."
  },
  {
    question: "Which toxin produced by Clostridium botulinum inhibits acetylcholine release at neuromuscular junctions, causing flaccid paralysis?",
    options: [
      "Tetanospasmin",
      "Botulinum toxin",
      "C. difficile toxin A",
      "Shiga toxin"
    ],
    answer: 1,
    explanation: "Botulinum toxin cleaves SNARE proteins in presynaptic terminals, preventing acetylcholine release and causing flaccid paralysis. Tetanospasmin (from C. tetani) causes spastic paralysis. C. difficile toxins cause colitis, and Shiga toxin inhibits protein synthesis in host cells."
  },
  {
    question: "Which bacterial pathogen uses actin-based motility to spread directly from cell to cell, avoiding extracellular exposure?",
    options: [
      "Listeria monocytogenes",
      "Salmonella enterica",
      "Streptococcus pyogenes",
      "Escherichia coli"
    ],
    answer: 0,
    explanation: "Listeria monocytogenes polymerizes host actin to propel itself into neighboring cells via membrane protrusions, enabling intracellular spread without extracellular exposure. Salmonella invades cells but does not use actin comet tails, Streptococcus and E. coli do not use this mechanism."
  },
  {
    question: "Which gram-negative pathogen causes dysentery via intracellular invasion and actin polymerization?",
    options: [
      "Shigella spp.",
      "Salmonella typhi",
      "Escherichia coli O157:H7",
      "Vibrio cholerae"
    ],
    answer: 0,
    explanation: "Shigella species invade colonic epithelial cells and use actin-based motility to move between cells, causing mucosal destruction and dysentery. Salmonella may invade intestines but causes systemic typhoid, EHEC produces Shiga toxin, and Vibrio cholerae causes watery diarrhea via cholera toxin."
  },
  {
    question: "Which factor contributes most to Pseudomonas aeruginosa's persistence in chronic infections?",
    options: [
      "Capsule polysaccharide alone",
      "Exotoxin A exclusively",
      "Biofilm formation and multidrug resistance",
      "Endotoxin only"
    ],
    answer: 2,
    explanation: "Biofilms protect P. aeruginosa from antibiotics and host immune responses, while multidrug resistance further complicates treatment. Exotoxin A and endotoxin contribute to virulence but are not sufficient for persistence."
  },
  {
    question: "Which bacterial virulence factor is most associated with Streptococcus pyogenes post-streptococcal rheumatic fever?",
    options: [
      "M protein",
      "Capsule",
      "Hyaluronidase",
      "Exotoxin A"
    ],
    answer: 0,
    explanation: "M protein exhibits molecular mimicry with cardiac tissue, leading to autoimmune sequelae like rheumatic fever. Capsule helps evade phagocytosis, hyaluronidase degrades connective tissue, and Exotoxin A acts as a superantigen."
  },
  {
    question: "Which molecular mechanism mediates antibiotic resistance in beta-lactamase-producing bacteria?",
    options: [
      "Efflux pumps only",
      "Enzymatic hydrolysis of beta-lactam ring",
      "Alteration of 16S rRNA",
      "Biofilm formation only"
    ],
    answer: 1,
    explanation: "Beta-lactamases hydrolyze the beta-lactam ring of penicillins and cephalosporins, inactivating the antibiotic. Efflux pumps, rRNA modification, and biofilms contribute to resistance but are distinct mechanisms."
  },
  {
    question: "Which viral protein of HIV is responsible for binding CD4 receptors and facilitating viral entry?",
    options: [
      "gp120",
      "p24",
      "Rev",
      "Tat"
    ],
    answer: 0,
    explanation: "gp120 on the HIV envelope binds CD4 and co-receptors (CCR5/CXCR4), initiating viral entry. p24 is a capsid protein; Rev and Tat are regulatory proteins involved in viral replication."
  },
  {
    question: "Which fungal pathogen exhibits dimorphism, existing as mold in the environment and yeast in host tissues?",
    options: [
      "Candida albicans",
      "Aspergillus fumigatus",
      "Histoplasma capsulatum",
      "Cryptococcus neoformans"
    ],
    answer: 2,
    explanation: "Histoplasma capsulatum grows as mold in soil and converts to yeast form in host at 37°C. Candida can form hyphae but is not strictly dimorphic, Aspergillus remains filamentous, Cryptococcus is encapsulated yeast."
  },
  {
    question: "Which protozoan pathogen is responsible for malaria and invades red blood cells?",
    options: [
      "Trypanosoma brucei",
      "Plasmodium falciparum",
      "Giardia lamblia",
      "Entamoeba histolytica"
    ],
    answer: 1,
    explanation: "Plasmodium falciparum is an intraerythrocytic parasite causing severe malaria. Trypanosoma brucei causes African sleeping sickness, Giardia affects intestines, Entamoeba causes dysentery."
  },
  {
    question: "Which helminth infection commonly modulates host immunity to evade clearance and persist chronically?",
    options: [
      "Enterobius vermicularis",
      "Schistosoma spp.",
      "Taenia solium",
      "Ascaris lumbricoides"
    ],
    answer: 1,
    explanation: "Schistosomes release immunomodulatory molecules that skew host response toward Th2 and regulatory pathways, promoting chronic infection. Other helminths may elicit immune responses but not as potent in immunomodulation."
  },
  {
    question: "Which exotoxin is responsible for the massive watery diarrhea seen in cholera?",
    options: [
      "Cholera toxin (CTX)",
      "Shiga toxin",
      "Botulinum toxin",
      "Tetanospasmin"
    ],
    answer: 0,
    explanation: "Cholera toxin is an A-B toxin that activates adenylate cyclase in intestinal epithelial cells, increasing cAMP and causing chloride and water secretion. Shiga toxin inhibits protein synthesis; botulinum and tetanus toxins affect neuromuscular function."
  },
  {
    question: "Which bacterial factor allows survival in acidic environments, such as the stomach?",
    options: [
      "Urease production",
      "Catalase",
      "Coagulase",
      "Protease"
    ],
    answer: 0,
    explanation: "Urease converts urea to ammonia and CO2, neutralizing gastric acid, enabling survival of Helicobacter pylori. Catalase breaks down hydrogen peroxide, coagulase promotes clotting, proteases degrade proteins but do not neutralize acid."
  },
  {
    question: "Which mechanism allows Listeria monocytogenes to escape the phagosome and replicate in the cytoplasm?",
    options: [
      "Listeriolysin O pore formation",
      "Endotoxin release",
      "Capsule formation",
      "Exotoxin A"
    ],
    answer: 0,
    explanation: "Listeriolysin O forms pores in the phagosomal membrane, allowing Listeria to enter the cytoplasm, avoiding lysosomal degradation. Endotoxin is a Gram-negative LPS component; capsule prevents phagocytosis, Exotoxin A is from Pseudomonas."
  },
  {
    question: "Which vector is primarily responsible for transmission of Plasmodium falciparum?",
    options: [
      "Aedes aegypti",
      "Anopheles mosquito",
      "Culex mosquito",
      "Tsetse fly"
    ],
    answer: 1,
    explanation: "Anopheles female mosquitoes transmit Plasmodium spp. to humans during blood feeding. Aedes transmits dengue/Zika, Culex transmits West Nile, Tsetse fly transmits Trypanosoma brucei."
  },
  {
    question: "Which molecular mechanism allows antigenic variation in Trypanosoma brucei?",
    options: [
      "VSG gene switching",
      "Phase variation of fimbriae",
      "Capsule synthesis",
      "Endospore formation"
    ],
    answer: 0,
    explanation: "Trypanosomes switch variant surface glycoprotein (VSG) genes periodically, evading host antibody responses. Phase variation applies to bacteria, capsules help immune evasion but not variation, endospores are bacterial survival structures."
  },
  {
    question: "Which pathogen's biofilm formation on medical devices contributes to chronic infections and antibiotic tolerance?",
    options: [
      "Staphylococcus epidermidis",
      "Clostridium tetani",
      "Neisseria gonorrhoeae",
      "Treponema pallidum"
    ],
    answer: 0,
    explanation: "S. epidermidis forms biofilms on catheters and prosthetic devices, protecting it from immune attack and antibiotics. Clostridium tetani produces tetanus toxin, Neisseria and Treponema rarely form device-associated biofilms."
  },
  {
    question: "Which bacterial pathogen produces Shiga toxin leading to hemolytic uremic syndrome (HUS)?",
    options: [
      "Escherichia coli O157:H7",
      "Salmonella typhi",
      "Vibrio cholerae",
      "Listeria monocytogenes"
    ],
    answer: 0,
    explanation: "EHEC strains of E. coli produce Shiga toxin that damages renal endothelium, leading to HUS. Salmonella causes typhoid fever, Vibrio cholerae causes watery diarrhea, Listeria causes systemic infection."
  },
  {
    question: "Which bacterial pathogen is associated with post-antibiotic colitis due to toxin production?",
    options: [
      "Clostridium difficile",
      "Streptococcus pyogenes",
      "Listeria monocytogenes",
      "Staphylococcus aureus"
    ],
    answer: 0,
    explanation: "C. difficile produces toxins A and B after dysbiosis from antibiotics, leading to colitis. Other bacteria do not typically cause post-antibiotic colitis."
  },
  {
    question: "Which fungal pathogen causes cryptococcal meningitis primarily in immunocompromised patients?",
    options: [
      "Candida albicans",
      "Cryptococcus neoformans",
      "Aspergillus fumigatus",
      "Histoplasma capsulatum"
    ],
    answer: 1,
    explanation: "Cryptococcus neoformans has a polysaccharide capsule that protects it from phagocytosis and allows dissemination to the CNS in immunocompromised patients. Candida and Aspergillus can cause systemic infections but less CNS involvement; Histoplasma mainly causes pulmonary disease."
  },
  {
    question: "Which bacterial pathogen is transmitted via ingestion of contaminated food or water and produces enterotoxin causing profuse watery diarrhea?",
    options: [
      "Vibrio cholerae",
      "Shigella dysenteriae",
      "Salmonella enterica",
      "Escherichia coli O157:H7"
    ],
    answer: 0,
    explanation: "Vibrio cholerae produces cholera toxin that activates adenylate cyclase in intestinal epithelial cells, causing massive fluid loss. Shigella and EHEC cause bloody diarrhea, Salmonella can cause gastroenteritis and systemic infection."
  },
  {
    question: "Which pathogen produces neurotoxin that prevents release of inhibitory neurotransmitters, causing spastic paralysis?",
    options: [
      "Clostridium tetani",
      "Clostridium botulinum",
      "Listeria monocytogenes",
      "Bacillus cereus"
    ],
    answer: 0,
    explanation: "Tetanospasmin from C. tetani inhibits release of GABA and glycine from inhibitory neurons, leading to uncontrolled muscle contractions (spastic paralysis). Botulinum toxin causes flaccid paralysis by inhibiting acetylcholine release."
  },
  {
    question: "Which bacterial virulence factor allows Helicobacter pylori to colonize the gastric mucosa?",
    options: [
      "Urease",
      "Capsule",
      "Type III secretion",
      "Lipoteichoic acid"
    ],
    answer: 0,
    explanation: "Urease converts urea to ammonia, neutralizing stomach acid and allowing H. pylori to survive in the acidic gastric environment. Capsule aids immune evasion, type III secretion is used by Gram-negative enteric pathogens, lipoteichoic acid is in Gram-positive bacteria."
  }
];


let quizData7 = [
  {
    question: "Which viral genome type requires RNA-dependent RNA polymerase for replication?",
    options: [
      "dsDNA",
      "(+)ssRNA",
      "(-)ssRNA",
      "Retrovirus (ssRNA-RT)"
    ],
    answer: 2,
    explanation: "(-)ssRNA viruses cannot be directly translated; they require viral RNA-dependent RNA polymerase to synthesize complementary (+)RNA for translation. (+)ssRNA can serve directly as mRNA, dsDNA uses host polymerases, retroviruses reverse transcribe RNA into DNA."
  },
  {
    question: "Which viral protein mediates attachment of HIV to CD4 and co-receptors on host cells?",
    options: ["gp120", "gp41", "p24", "Tat"],
    answer: 0,
    explanation: "gp120 binds CD4 and co-receptors CCR5/CXCR4, initiating viral entry. gp41 mediates fusion, p24 is capsid, Tat is a transcriptional regulator."
  },
  {
    question: "Which virus family is known for latency in neurons and reactivation causing cold sores?",
    options: ["Herpesviridae", "Poxviridae", "Adenoviridae", "Picornaviridae"],
    answer: 0,
    explanation: "Herpes simplex virus (HSV) establishes latency in sensory neurons and can reactivate causing cold sores. Poxviruses, adenoviruses, and picornaviruses do not establish neuronal latency."
  },
  {
    question: "Which viral enzyme is targeted by nucleoside analogues like acyclovir?",
    options: ["Reverse transcriptase", "DNA polymerase", "RNA polymerase", "Integrase"],
    answer: 1,
    explanation: "Acyclovir is phosphorylated by viral thymidine kinase and inhibits viral DNA polymerase, preventing DNA synthesis. Reverse transcriptase is targeted by HIV drugs, RNA polymerase by RNA virus antivirals, integrase by HIV integrase inhibitors."
  },
  {
    question: "Which virus has a segmented RNA genome allowing antigenic shift?",
    options: ["Influenza A virus", "Hepatitis B virus", "Dengue virus", "HIV"],
    answer: 0,
    explanation: "Influenza A has 8 segmented (-)ssRNA segments; reassortment between strains causes antigenic shift. Hepatitis B is DNA, Dengue is (+)ssRNA non-segmented, HIV is retrovirus."
  },
  {
    question: "Which viral family includes enveloped (+)ssRNA viruses responsible for Zika and Dengue?",
    options: ["Flaviviridae", "Picornaviridae", "Retroviridae", "Orthomyxoviridae"],
    answer: 0,
    explanation: "Flaviviridae contains enveloped (+)ssRNA viruses like Dengue, Zika, and Yellow fever. Picornaviridae are non-enveloped, Retroviridae are ssRNA-RT, Orthomyxoviridae are segmented (-)ssRNA."
  },
  {
    question: "Which viral diagnostic method detects viral nucleic acid quantitatively?",
    options: ["ELISA", "qPCR", "Hemagglutination assay", "Culture"],
    answer: 1,
    explanation: "qPCR quantifies viral nucleic acids using fluorescent probes during amplification. ELISA detects antibodies/antigens, hemagglutination detects viral particles functionally, culture shows cytopathic effects."
  },
  {
    question: "Which viral strategy allows evasion of neutralizing antibodies via frequent mutations in surface glycoproteins?",
    options: ["Antigenic variation", "Latency", "Lytic replication", "Integration"],
    answer: 0,
    explanation: "Antigenic variation enables viruses (like influenza) to evade immune recognition by altering epitopes. Latency hides genome without protein expression, lytic replication kills host cells, integration inserts viral DNA into host genome."
  },
  {
    question: "Which virus uses reverse transcription to integrate into the host genome?",
    options: ["HIV", "Influenza A", "Hepatitis C", "Herpes simplex virus"],
    answer: 0,
    explanation: "HIV (Retroviridae) converts its ssRNA genome into DNA using reverse transcriptase, then integrates into host chromosomes. Influenza, HCV, and HSV do not undergo reverse transcription."
  },
  {
    question: "Which viral protein mediates fusion of HIV envelope with host cell membrane?",
    options: ["gp41", "gp120", "p24", "Rev"],
    answer: 0,
    explanation: "gp41 facilitates membrane fusion after gp120 binds CD4/co-receptor. gp120 binds receptor, p24 forms capsid, Rev regulates transcription."
  },
  {
    question: "Which virus causes chronic liver infection and is a DNA virus with a reverse transcription step?",
    options: ["Hepatitis B virus", "Hepatitis C virus", "Hepatitis A virus", "Hepatitis E virus"],
    answer: 0,
    explanation: "Hepatitis B virus is a partially double-stranded DNA virus that uses reverse transcription during replication. HCV is RNA, HAV/HEV are RNA viruses without DNA stage."
  },
  {
    question: "Which viral replication strategy occurs entirely in the host cell cytoplasm?",
    options: ["Poxvirus", "Herpesvirus", "Adenovirus", "Papillomavirus"],
    answer: 0,
    explanation: "Poxviruses replicate entirely in cytoplasm because they carry all necessary enzymes for transcription and replication. Herpes, adenoviruses, papillomaviruses replicate in the nucleus."
  },
  {
    question: "Which viral diagnostic test measures host antibody response to infection?",
    options: ["PCR", "qPCR", "ELISA", "Electron microscopy"],
    answer: 2,
    explanation: "ELISA detects antibodies produced by the host against viral antigens, indicating past or current infection. PCR/qPCR detect viral nucleic acids directly, EM visualizes virus particles."
  },
  {
    question: "Which viral component is most important for determining host cell tropism?",
    options: ["Capsid proteins", "Envelope glycoproteins", "Polymerase", "Integrase"],
    answer: 1,
    explanation: "Envelope glycoproteins bind specific host receptors, determining which cells the virus can infect. Capsid protects genome, polymerase is for replication, integrase inserts viral DNA."
  },
  {
    question: "Which viral family is non-enveloped and causes common cold?",
    options: ["Picornaviridae", "Flaviviridae", "Retroviridae", "Herpesviridae"],
    answer: 0,
    explanation: "Picornaviridae (e.g., Rhinovirus) are non-enveloped and cause mild respiratory infections. Flaviviridae are enveloped, Retroviridae and Herpesviridae are enveloped DNA/RNA viruses."
  },
  {
    question: "Which antiviral drug class targets influenza virus neuraminidase?",
    options: ["Oseltamivir", "Acyclovir", "Zidovudine", "Ribavirin"],
    answer: 0,
    explanation: "Oseltamivir inhibits neuraminidase, preventing release of progeny influenza virus. Acyclovir targets viral DNA polymerase, Zidovudine targets HIV reverse transcriptase, Ribavirin inhibits viral RNA synthesis broadly."
  },
  {
    question: "Which virus is an example of an enveloped, negative-sense, single-stranded RNA virus with segmented genome?",
    options: ["Influenza A", "Hepatitis B", "Dengue virus", "HIV"],
    answer: 0,
    explanation: "Influenza A virus has (-)ssRNA segmented genome and envelope. HBV is DNA, Dengue is (+)ssRNA enveloped, HIV is ssRNA-RT."
  },
  {
    question: "Which antiviral strategy uses monoclonal antibodies to neutralize virus particles?",
    options: ["Passive immunization", "Active immunization", "Nucleoside analogues", "Protease inhibitors"],
    answer: 0,
    explanation: "Passive immunization delivers preformed antibodies to neutralize virus immediately. Active immunization stimulates host immunity, nucleoside analogues inhibit replication, protease inhibitors block viral enzymes."
  },
  {
    question: "Which viral disease is primarily controlled through vector management rather than vaccination?",
    options: ["Dengue", "Measles", "Hepatitis B", "Polio"],
    answer: 0,
    explanation: "Dengue is transmitted by Aedes mosquitoes; vector control is critical. Measles, Hepatitis B, and Polio are controlled mainly through vaccination."
  },
  {
    question: "Which virus uses antigenic shift as a mechanism for emergence of pandemic strains?",
    options: ["Influenza A", "Hepatitis C", "HIV", "Rabies virus"],
    answer: 0,
    explanation: "Influenza A undergoes reassortment of segmented RNA genomes between strains (antigenic shift), potentially causing pandemics. HCV mutates slowly, HIV uses antigenic variation, rabies virus is stable."
  },
  {
    question: "Which viral infection is characterized by lifelong latent infection in neurons with periodic reactivation?",
    options: ["Herpes simplex virus", "Hepatitis A virus", "Poliovirus", "Adenovirus"],
    answer: 0,
    explanation: "Herpes simplex virus establishes latency in sensory neurons and can reactivate to cause cold sores or genital lesions. Hepatitis A, Poliovirus, Adenovirus do not establish neuronal latency."
  },
  {
    question: "Which virus family contains single-stranded RNA viruses that replicate via a DNA intermediate?",
    options: ["Retroviridae", "Flaviviridae", "Picornaviridae", "Orthomyxoviridae"],
    answer: 0,
    explanation: "Retroviridae (e.g., HIV) use reverse transcription of ssRNA into DNA for integration into host genome. Flaviviridae, Picornaviridae, and Orthomyxoviridae replicate RNA without DNA intermediate."
  },
  {
    question: "Which virus can cause severe congenital defects if infection occurs during pregnancy?",
    options: ["Zika virus", "Rhinovirus", "Influenza A", "Adenovirus"],
    answer: 0,
    explanation: "Zika virus infection during pregnancy can cause microcephaly and other congenital malformations. Rhinovirus causes mild cold, Influenza A rarely causes congenital defects, Adenovirus mainly respiratory disease."
  },
  {
    question: "Which method allows direct visualization of viral particles for identification?",
    options: ["Electron microscopy", "PCR", "ELISA", "qPCR"],
    answer: 0,
    explanation: "Electron microscopy provides direct visualization of virus morphology. PCR/qPCR detect nucleic acids, ELISA detects antibodies or antigens."
  }
];

let quizData8 = [
  {
    question: "Which component of the fungal cell membrane is targeted by azole antifungals?",
    options: ["Chitin", "Ergosterol", "Glucan", "Mannoproteins"],
    answer: 1,
    explanation: "Azoles inhibit lanosterol 14α-demethylase, blocking ergosterol synthesis, which disrupts the fungal cell membrane. Chitin is part of the cell wall, glucan is cell wall polysaccharide, mannoproteins are surface proteins."
  },
  {
    question: "Which dimorphic fungus grows as mold at 25°C and yeast at 37°C, causing histoplasmosis?",
    options: ["Candida albicans", "Histoplasma capsulatum", "Aspergillus fumigatus", "Cryptococcus neoformans"],
    answer: 1,
    explanation: "Histoplasma capsulatum is a dimorphic fungus: mold in the environment (25°C) and yeast in host tissue (37°C). Candida is mostly yeast, Aspergillus is filamentous, Cryptococcus is encapsulated yeast."
  },
  {
    question: "Which fungus produces a thick polysaccharide capsule allowing immune evasion and CNS infection?",
    options: ["Candida albicans", "Aspergillus fumigatus", "Cryptococcus neoformans", "Histoplasma capsulatum"],
    answer: 2,
    explanation: "Cryptococcus neoformans has a polysaccharide capsule that inhibits phagocytosis and complement activation, facilitating meningitis, especially in immunocompromised hosts."
  },
  {
    question: "Which fungal reproductive structure is asexual in Aspergillus spp.?",
    options: ["Zygospore", "Conidium", "Ascospore", "Basidiospore"],
    answer: 1,
    explanation: "Aspergillus produces conidia (asexual spores) that disperse in air. Zygospore is sexual (Zygomycota), ascospores (Ascomycota), basidiospores (Basidiomycota)."
  },
  {
    question: "Which antifungal class inhibits β-(1,3)-D-glucan synthesis in the fungal cell wall?",
    options: ["Polyenes", "Azoles", "Echinocandins", "Flucytosine"],
    answer: 2,
    explanation: "Echinocandins inhibit β-(1,3)-D-glucan synthase, compromising the fungal cell wall. Polyenes bind ergosterol, azoles inhibit ergosterol synthesis, flucytosine inhibits nucleic acid synthesis."
  },
  {
    question: "Which dermatophyte genus commonly infects hair, skin, and nails?",
    options: ["Trichophyton", "Aspergillus", "Candida", "Cryptococcus"],
    answer: 0,
    explanation: "Trichophyton species infect keratinized tissues, causing tinea infections. Aspergillus causes systemic or respiratory infections, Candida causes mucosal and systemic infections, Cryptococcus causes CNS infection."
  },
  {
    question: "Which mycotoxin is produced by Aspergillus flavus and is hepatotoxic and carcinogenic?",
    options: ["Ochratoxin", "Aflatoxin", "Zearalenone", "Fumonisin"],
    answer: 1,
    explanation: "Aflatoxin B1 produced by Aspergillus flavus is hepatotoxic and hepatocarcinogenic. Ochratoxin is nephrotoxic, Zearalenone is estrogenic, Fumonisin affects sphingolipid metabolism."
  },
  {
    question: "Which fungal diagnostic method uses KOH to dissolve human tissue and visualize fungal hyphae?",
    options: ["Culture", "Microscopy", "Serology", "PCR"],
    answer: 1,
    explanation: "KOH mount clears keratin and human tissue, allowing visualization of fungal elements under microscope. Culture grows fungus, serology detects antibodies/antigens, PCR detects nucleic acids."
  },
  {
    question: "Which fungal pathogen forms septate hyphae with acute-angle branching, causing invasive pulmonary infections?",
    options: ["Aspergillus spp.", "Mucor spp.", "Candida albicans", "Cryptococcus neoformans"],
    answer: 0,
    explanation: "Aspergillus produces septate hyphae with 45° (acute-angle) branching. Mucor has broad non-septate hyphae, Candida forms yeast or pseudohyphae, Cryptococcus is yeast."
  },
  {
    question: "Which opportunistic yeast is a major cause of systemic infections in immunocompromised patients?",
    options: ["Candida albicans", "Aspergillus fumigatus", "Trichophyton", "Histoplasma capsulatum"],
    answer: 0,
    explanation: "Candida albicans can cause mucosal and systemic infections (candidemia) in immunocompromised patients. Aspergillus causes respiratory or systemic infections, Trichophyton causes superficial infections, Histoplasma causes dimorphic systemic infections."
  },
  {
    question: "Which fungal structure is formed by sexual reproduction in Zygomycota?",
    options: ["Zygospore", "Conidium", "Ascospore", "Basidiospore"],
    answer: 0,
    explanation: "Zygomycota reproduce sexually via zygospores. Conidia are asexual spores, ascospores from Ascomycota, basidiospores from Basidiomycota."
  },
  {
    question: "Which fungal enzyme degrades keratin in hair and nails?",
    options: ["Protease", "Keratinase", "Lipase", "Amylase"],
    answer: 1,
    explanation: "Keratinase specifically degrades keratin, enabling dermatophytes to invade hair, skin, and nails. Protease is general, lipase digests fats, amylase digests starch."
  },
  {
    question: "Which antifungal drug inhibits nucleic acid synthesis in fungi by conversion to 5-fluorouracil?",
    options: ["Flucytosine", "Amphotericin B", "Fluconazole", "Caspofungin"],
    answer: 0,
    explanation: "Flucytosine enters fungal cells and is converted to 5-FU, inhibiting DNA/RNA synthesis. Amphotericin binds ergosterol, fluconazole inhibits ergosterol synthesis, caspofungin inhibits glucan synthesis."
  },
  {
    question: "Which fungal pathogen causes systemic infection via inhalation of spores and is thermally dimorphic?",
    options: ["Histoplasma capsulatum", "Candida albicans", "Aspergillus fumigatus", "Cryptococcus neoformans"],
    answer: 0,
    explanation: "Histoplasma capsulatum spores are inhaled, convert to yeast at 37°C, causing pulmonary and systemic histoplasmosis. Candida is opportunistic yeast, Aspergillus is filamentous, Cryptococcus is yeast."
  },
  {
    question: "Which fungal diagnostic stain is used to visualize capsules of Cryptococcus?",
    options: ["India ink", "PAS", "Calcofluor white", "Giemsa"],
    answer: 0,
    explanation: "India ink provides a negative stain that highlights the polysaccharide capsule of Cryptococcus. PAS stains fungal cell walls, Calcofluor binds chitin, Giemsa is used for blood parasites."
  },
  {
    question: "Which type of hyphae lack septa and are typical of Mucorales?",
    options: ["Septate hyphae", "Coenocytic hyphae", "Pseudohyphae", "Yeast cells"],
    answer: 1,
    explanation: "Coenocytic (non-septate) hyphae are characteristic of Mucorales. Septate hyphae are in Aspergillus, pseudohyphae in Candida, yeast are unicellular."
  },
  {
    question: "Which antifungal class binds ergosterol and forms pores causing fungal cell death?",
    options: ["Polyenes", "Azoles", "Echinocandins", "Flucytosine"],
    answer: 0,
    explanation: "Polyenes (Amphotericin B) bind ergosterol in fungal membranes, forming pores that disrupt integrity, leading to cell death."
  },
  {
    question: "Which superficial fungal infection affects the scalp and hair shafts?",
    options: ["Tinea capitis", "Candidiasis", "Histoplasmosis", "Aspergillosis"],
    answer: 0,
    explanation: "Tinea capitis is caused by dermatophytes (Trichophyton, Microsporum) infecting hair and scalp. Candidiasis affects mucosa or skin, Histoplasmosis is systemic, Aspergillosis mainly respiratory."
  },
  {
    question: "Which environmental condition favors fungal growth in stored food, increasing mycotoxin production?",
    options: ["High moisture and warmth", "Low moisture", "Freezing", "Dry air"],
    answer: 0,
    explanation: "High moisture and warm temperatures favor fungal growth and mycotoxin production. Low moisture, freezing, and dry air inhibit fungal proliferation."
  },
  {
    question: "Which fungal pathogen produces pseudohyphae and true hyphae, contributing to pathogenicity?",
    options: ["Candida albicans", "Aspergillus fumigatus", "Cryptococcus neoformans", "Trichophyton spp."],
    answer: 0,
    explanation: "Candida albicans exhibits dimorphism with yeast, pseudohyphae, and hyphae; hyphal forms invade tissues and enhance virulence. Aspergillus is filamentous, Cryptococcus is yeast, Trichophyton is septate mold."
  },
  {
    question: "Which immunocompromised population is most at risk for invasive Aspergillus infection?",
    options: ["Neutropenic patients", "Healthy adults", "Children with common cold", "Pregnant women"],
    answer: 0,
    explanation: "Neutropenic patients lack effective innate immunity, increasing susceptibility to invasive Aspergillus infections. Healthy adults and children rarely develop invasive disease."
  },
  {
    question: "Which fungal structure is used for sexual reproduction in Ascomycota?",
    options: ["Ascospore", "Basidiospore", "Zygospore", "Conidium"],
    answer: 0,
    explanation: "Ascomycota produce ascospores in asci during sexual reproduction. Basidiospores are Basidiomycota, zygospores are Zygomycota, conidia are asexual spores."
  },
  {
    question: "Which fungal pathogen commonly causes oropharyngeal and vaginal infections in humans?",
    options: ["Candida albicans", "Aspergillus fumigatus", "Cryptococcus neoformans", "Histoplasma capsulatum"],
    answer: 0,
    explanation: "Candida albicans causes mucosal infections such as thrush and vaginal candidiasis. Aspergillus and Histoplasma cause systemic or respiratory infections; Cryptococcus causes CNS infection."
  },
  {
    question: "Which laboratory method uses Sabouraud dextrose agar to culture fungi?",
    options: ["Fungal culture", "KOH mount", "PCR", "ELISA"],
    answer: 0,
    explanation: "Sabouraud dextrose agar supports growth of fungi in culture. KOH mount is for direct microscopy, PCR detects DNA, ELISA detects antigens or antibodies."
  }
];

let quizData9 = [
  {
    question: "Which stage of Plasmodium is responsible for liver infection in humans?",
    options: ["Sporozoite", "Merozoite", "Gametocyte", "Trophozoite"],
    answer: 0,
    explanation: "Sporozoites injected by Anopheles mosquitoes invade hepatocytes, initiating the liver stage. Merozoites infect red blood cells, gametocytes are sexual forms, trophozoites are intraerythrocytic forms."
  },
  {
    question: "Which protozoan causes watery diarrhea and has cyst and trophozoite stages?",
    options: ["Giardia lamblia", "Plasmodium falciparum", "Entamoeba histolytica", "Trypanosoma brucei"],
    answer: 0,
    explanation: "Giardia lamblia causes giardiasis with cysts (infective) and trophozoites (pathogenic) stages. Plasmodium infects RBCs, Entamoeba causes dysentery, Trypanosoma causes sleeping sickness."
  },
  {
    question: "Which helminth is transmitted via ingestion of eggs in contaminated soil and causes intestinal obstruction?",
    options: ["Ascaris lumbricoides", "Schistosoma mansoni", "Enterobius vermicularis", "Ancylostoma duodenale"],
    answer: 0,
    explanation: "Ascaris lumbricoides eggs are ingested, larvae migrate through tissues and mature in the intestines, potentially causing obstruction. Schistosoma penetrates skin, Enterobius causes perianal itching, Ancylostoma penetrates skin causing anemia."
  },
  {
    question: "Which parasite forms cysts in the brain causing toxoplasmosis in immunocompromised patients?",
    options: ["Toxoplasma gondii", "Plasmodium falciparum", "Giardia lamblia", "Entamoeba histolytica"],
    answer: 0,
    explanation: "Toxoplasma gondii forms tissue cysts in the brain and muscles; immunocompromised hosts may develop encephalitis. Plasmodium affects RBCs, Giardia causes intestinal infection, Entamoeba invades colon."
  },
  {
    question: "Which stage of Plasmodium is taken up by the mosquito to continue the life cycle?",
    options: ["Gametocyte", "Sporozoite", "Merozoite", "Trophozoite"],
    answer: 0,
    explanation: "Gametocytes are sexual stages ingested by Anopheles mosquitoes to continue the parasite's life cycle. Sporozoites infect humans, merozoites infect RBCs, trophozoites are intraerythrocytic stages."
  },
  {
    question: "Which helminth infection is characterized by perianal itching and is diagnosed using the Scotch tape test?",
    options: ["Enterobius vermicularis", "Ascaris lumbricoides", "Trichuris trichiura", "Ancylostoma duodenale"],
    answer: 0,
    explanation: "Enterobius vermicularis (pinworm) lays eggs around the perianal region; Scotch tape test collects eggs for microscopy. Ascaris and Trichuris eggs are in stool; Ancylostoma penetrates skin."
  },
  {
    question: "Which parasite is transmitted via sandfly bite and causes visceral or cutaneous leishmaniasis?",
    options: ["Leishmania spp.", "Plasmodium spp.", "Trypanosoma brucei", "Giardia lamblia"],
    answer: 0,
    explanation: "Leishmania promastigotes are transmitted by sandflies, causing cutaneous or visceral disease. Plasmodium is mosquito-borne, Trypanosoma brucei by tsetse fly, Giardia via fecal-oral route."
  },
  {
    question: "Which trematode requires a snail as an intermediate host and causes urinary or intestinal schistosomiasis?",
    options: ["Schistosoma spp.", "Fasciola hepatica", "Clonorchis sinensis", "Paragonimus westermani"],
    answer: 0,
    explanation: "Schistosoma spp. use freshwater snails as intermediate hosts. Fasciola infects liver via ingestion, Clonorchis infects bile ducts, Paragonimus infects lungs."
  },
  {
    question: "Which parasite causes sleeping sickness and is transmitted by tsetse flies?",
    options: ["Trypanosoma brucei", "Plasmodium falciparum", "Leishmania donovani", "Giardia lamblia"],
    answer: 0,
    explanation: "Trypanosoma brucei is transmitted by Glossina tsetse flies, causing African trypanosomiasis (sleeping sickness). Plasmodium causes malaria, Leishmania causes leishmaniasis, Giardia causes giardiasis."
  },
  {
    question: "Which intestinal protozoan invades colonic mucosa causing dysentery and flask-shaped ulcers?",
    options: ["Entamoeba histolytica", "Giardia lamblia", "Plasmodium falciparum", "Toxoplasma gondii"],
    answer: 0,
    explanation: "Entamoeba histolytica trophozoites invade colonic mucosa causing flask-shaped ulcers and bloody diarrhea. Giardia causes watery diarrhea, Plasmodium infects RBCs, Toxoplasma infects tissues."
  },
  {
    question: "Which parasite causes lymphatic filariasis and is transmitted by mosquitoes?",
    options: ["Wuchereria bancrofti", "Ascaris lumbricoides", "Enterobius vermicularis", "Schistosoma mansoni"],
    answer: 0,
    explanation: "Wuchereria bancrofti microfilariae are transmitted by mosquitoes, causing lymphatic obstruction and elephantiasis. Ascaris and Enterobius are intestinal worms, Schistosoma causes schistosomiasis."
  },
  {
    question: "Which protozoan undergoes antigenic variation to evade host immunity in chronic infections?",
    options: ["Trypanosoma brucei", "Giardia lamblia", "Plasmodium falciparum", "Entamoeba histolytica"],
    answer: 0,
    explanation: "Trypanosoma brucei expresses variable surface glycoproteins to evade host immunity. Plasmodium changes antigens but less dynamically, Giardia and Entamoeba do not have high antigenic variation."
  },
  {
    question: "Which cestode infection can lead to neurocysticercosis if humans ingest eggs?",
    options: ["Taenia solium", "Taenia saginata", "Diphyllobothrium latum", "Echinococcus granulosus"],
    answer: 0,
    explanation: "Taenia solium eggs ingested by humans can develop into cysticerci in CNS, causing neurocysticercosis. T. saginata mainly causes intestinal infection, Diphyllobothrium affects vitamin B12 absorption, Echinococcus causes hydatid cysts."
  },
  {
    question: "Which parasite is diagnosed by thick and thin blood smears showing ring forms in RBCs?",
    options: ["Plasmodium spp.", "Trypanosoma brucei", "Leishmania spp.", "Giardia lamblia"],
    answer: 0,
    explanation: "Plasmodium trophozoites appear as ring forms in RBCs; thick smears detect parasites, thin smears allow species identification. Trypanosoma is in blood as trypomastigotes, Leishmania in tissue, Giardia in stool."
  },
  {
    question: "Which helminth infection causes iron-deficiency anemia via blood feeding in the intestine?",
    options: ["Ancylostoma duodenale", "Ascaris lumbricoides", "Enterobius vermicularis", "Trichuris trichiura"],
    answer: 0,
    explanation: "Hookworms (Ancylostoma duodenale, Necator americanus) attach to intestinal mucosa and feed on blood, causing anemia. Ascaris may cause obstruction, Enterobius causes itching, Trichuris causes diarrhea."
  },
  {
    question: "Which protozoan is transmitted via cysts in contaminated water and causes malabsorption?",
    options: ["Giardia lamblia", "Entamoeba histolytica", "Plasmodium falciparum", "Toxoplasma gondii"],
    answer: 0,
    explanation: "Giardia lamblia cysts in contaminated water cause malabsorption and chronic diarrhea. Entamoeba invades colon, Plasmodium infects RBCs, Toxoplasma forms tissue cysts."
  },
  {
    question: "Which trematode infects humans via ingestion of metacercariae in undercooked fish?",
    options: ["Clonorchis sinensis", "Schistosoma mansoni", "Fasciola hepatica", "Paragonimus westermani"],
    answer: 0,
    explanation: "Clonorchis sinensis metacercariae in freshwater fish are ingested, causing hepatobiliary disease. Schistosoma penetrates skin, Fasciola infects liver via watercress, Paragonimus infects lungs via crab/crayfish."
  },
  {
    question: "Which parasite causes visceral leishmaniasis, characterized by hepatosplenomegaly and pancytopenia?",
    options: ["Leishmania donovani", "Trypanosoma brucei", "Giardia lamblia", "Entamoeba histolytica"],
    answer: 0,
    explanation: "Leishmania donovani infects macrophages in liver, spleen, and bone marrow causing kala-azar (visceral leishmaniasis). Trypanosoma causes sleeping sickness, Giardia causes intestinal disease, Entamoeba causes colitis."
  },
  {
    question: "Which parasitic infection is commonly diagnosed using serological detection of specific IgG antibodies?",
    options: ["Toxoplasma gondii", "Plasmodium falciparum", "Giardia lamblia", "Enterobius vermicularis"],
    answer: 0,
    explanation: "Toxoplasmosis is often diagnosed by detecting IgG antibodies indicating past exposure or IgM for acute infection. Plasmodium is mainly diagnosed by blood smear, Giardia and Enterobius by stool."
  },
  {
    question: "Which intestinal nematode has a direct life cycle without intermediate hosts and causes whipworm infection?",
    options: ["Trichuris trichiura", "Ascaris lumbricoides", "Ancylostoma duodenale", "Enterobius vermicularis"],
    answer: 0,
    explanation: "Trichuris trichiura eggs are ingested from contaminated soil, hatch in intestine, causing whipworm infection. Ascaris also has direct cycle but causes roundworm infection, Ancylostoma penetrates skin, Enterobius causes pinworm."
  },
  {
    question: "Which parasite can cross the placenta and infect the fetus, leading to congenital infection?",
    options: ["Toxoplasma gondii", "Plasmodium falciparum", "Giardia lamblia", "Schistosoma mansoni"],
    answer: 0,
    explanation: "Toxoplasma gondii tachyzoites can cross the placenta causing congenital toxoplasmosis. Plasmodium can cause maternal malaria complications, Giardia and Schistosoma rarely cross placenta."
  },
  {
    question: "Which protozoan is responsible for American trypanosomiasis (Chagas disease) transmitted by triatomine bugs?",
    options: ["Trypanosoma cruzi", "Trypanosoma brucei", "Leishmania spp.", "Giardia lamblia"],
    answer: 0,
    explanation: "Trypanosoma cruzi is transmitted by triatomine bugs and causes Chagas disease affecting heart and gastrointestinal system. Trypanosoma brucei causes African sleeping sickness, Leishmania causes leishmaniasis, Giardia causes intestinal infection."
  },
  {
    question: "Which parasite's eggs are detected in stool using concentration techniques like flotation or sedimentation?",
    options: ["Helminths", "Protozoa", "Arthropods", "Viruses"],
    answer: 0,
    explanation: "Helminth eggs are detected in stool using flotation or sedimentation concentration methods for diagnosis. Protozoa are detected as cysts/trophozoites, arthropods and viruses are not detected in stool in this way."
  },
  {
    question: "Which parasitic infection can be prevented by boiling or filtering drinking water?",
    options: ["Giardiasis", "Malaria", "Leishmaniasis", "Filariasis"],
    answer: 0,
    explanation: "Giardia lamblia cysts in contaminated water can be killed by boiling or removed by filtration. Malaria and Leishmaniasis are vector-borne, Filariasis is mosquito-borne."
  }
];


let quizData10 = [
  {
    question: "Which term describes the constant presence of a disease within a population?",
    options: ["Epidemic", "Endemic", "Pandemic", "Sporadic"],
    answer: 1,
    explanation: "An endemic disease is constantly present in a population, such as malaria in certain tropical regions. Epidemic is a sudden increase, pandemic is global spread, sporadic occurs irregularly."
  },
  {
    question: "Which measure represents the number of new cases in a population over a specified period?",
    options: ["Prevalence", "Incidence", "Morbidity", "Mortality"],
    answer: 1,
    explanation: "Incidence refers to new cases over a time period, whereas prevalence includes both new and existing cases."
  },
  {
    question: "Which is the primary vector for transmitting Plasmodium spp.?",
    options: ["Aedes mosquito", "Anopheles mosquito", "Culex mosquito", "Tsetse fly"],
    answer: 1,
    explanation: "Anopheles mosquitoes transmit Plasmodium spp. causing malaria. Aedes transmits dengue and Zika, Culex transmits West Nile virus, Tsetse fly transmits Trypanosoma brucei."
  },
  {
    question: "Which term defines an individual who harbors a pathogen without showing symptoms but can transmit it?",
    options: ["Reservoir", "Carrier", "Vector", "Host"],
    answer: 1,
    explanation: "A carrier harbors a pathogen asymptomatically but can transmit it. A reservoir is any habitat, a vector transmits pathogens, a host is affected by the pathogen."
  },
  {
    question: "Which type of surveillance involves proactive collection of data by health authorities?",
    options: ["Passive surveillance", "Active surveillance", "Sentinel surveillance", "Syndromic surveillance"],
    answer: 1,
    explanation: "Active surveillance involves proactive efforts by public health authorities to identify cases. Passive relies on routine reporting, sentinel uses selected sites, syndromic uses symptom patterns."
  },
  {
    question: "Which intervention reduces vector-borne diseases by killing mosquitoes or reducing breeding sites?",
    options: ["Vaccination", "Vector control", "Sanitation", "Quarantine"],
    answer: 1,
    explanation: "Vector control targets insects or arthropods that transmit pathogens, including insecticides and removing breeding sites. Vaccination prevents infection, sanitation reduces fecal-oral spread, quarantine isolates infected individuals."
  },
  {
    question: "Which epidemiologic study design compares exposure in diseased and non-diseased groups to identify risk factors?",
    options: ["Cohort study", "Case-control study", "Cross-sectional study", "Randomized trial"],
    answer: 1,
    explanation: "Case-control studies compare individuals with disease (cases) and without (controls) to identify past exposures. Cohort follows exposed and unexposed, cross-sectional measures prevalence, randomized trial tests interventions."
  },
  {
    question: "Which step in outbreak investigation involves defining what constitutes a case?",
    options: ["Hypothesis testing", "Case definition", "Control measures", "Data collection"],
    answer: 1,
    explanation: "Defining a case is crucial to identify affected individuals accurately. Hypothesis testing and control measures come later, data collection supports both."
  },
  {
    question: "Which microorganism is commonly used as an indicator of water contamination in public health microbiology?",
    options: ["Escherichia coli", "Salmonella typhi", "Vibrio cholerae", "Clostridium difficile"],
    answer: 0,
    explanation: "E. coli presence indicates fecal contamination in water. Salmonella typhi and Vibrio cholerae are pathogens, not indicators, Clostridium difficile is healthcare-associated."
  },
  {
    question: "Which type of epidemiological curve helps identify the pattern and source of an outbreak?",
    options: ["Line listing", "Histogram", "Epidemic curve", "Scatter plot"],
    answer: 2,
    explanation: "An epidemic curve graphs cases over time, revealing patterns like point-source or propagated outbreaks."
  },
  {
    question: "Which term describes a sudden increase in disease cases above the expected level in a population?",
    options: ["Endemic", "Epidemic", "Pandemic", "Outbreak"],
    answer: 1,
    explanation: "An epidemic is a sudden rise above expected levels. Endemic is constant, pandemic is global, outbreak is localized epidemic."
  },
  {
    question: "Which strategy is essential for preventing foodborne outbreaks of Salmonella and E. coli?",
    options: ["Vector control", "Vaccination", "Food hygiene", "Quarantine"],
    answer: 2,
    explanation: "Proper food hygiene, cooking, and storage prevent contamination and outbreaks. Vector control targets insects, vaccination prevents infection, quarantine isolates cases."
  },
  {
    question: "Which measure is calculated as total cases divided by total population at a point in time?",
    options: ["Incidence", "Prevalence", "Attack rate", "Mortality rate"],
    answer: 1,
    explanation: "Prevalence measures all existing cases at a specific time, while incidence measures new cases. Attack rate is during outbreaks, mortality rate measures deaths."
  },
  {
    question: "Which pathogen is commonly monitored in hospital microbiology for nosocomial infection control?",
    options: ["Staphylococcus aureus", "Giardia lamblia", "Plasmodium falciparum", "Toxoplasma gondii"],
    answer: 0,
    explanation: "S. aureus, especially MRSA, is a common nosocomial pathogen. Giardia, Plasmodium, and Toxoplasma are mostly community-acquired."
  },
  {
    question: "Which type of disease transmission involves contaminated fomites or medical equipment?",
    options: ["Direct", "Indirect", "Vector-borne", "Airborne"],
    answer: 1,
    explanation: "Indirect transmission occurs via fomites or inanimate objects, unlike direct contact, vector-borne, or airborne transmission."
  },
  {
    question: "Which factor is considered part of the epidemiologic triad alongside host and agent?",
    options: ["Environment", "Vector", "Pathogenicity", "Transmission mode"],
    answer: 0,
    explanation: "The epidemiologic triad consists of host, agent, and environment which includes factors influencing disease occurrence. Vector is part of environment, pathogenicity is agent property, transmission mode is mechanism."
  },
  {
    question: "Which term describes diseases that reappear after being controlled due to environmental or social changes?",
    options: ["Emerging diseases", "Re-emerging diseases", "Endemic diseases", "Sporadic diseases"],
    answer: 1,
    explanation: "Re-emerging diseases reappear due to factors like reduced vaccination, resistance, or ecological changes. Emerging diseases are new infections, endemic are constant, sporadic are irregular."
  },
  {
    question: "Which microbiological test is commonly used to confirm cholera outbreaks in public health surveillance?",
    options: ["Stool culture for Vibrio cholerae", "Blood smear for Plasmodium", "KOH mount", "PCR for Toxoplasma"],
    answer: 0,
    explanation: "Stool culture detects Vibrio cholerae, confirming cholera outbreaks. Blood smear is for malaria, KOH mount for fungi, PCR for Toxoplasma DNA."
  },
  {
    question: "Which measure describes the proportion of exposed individuals who develop disease in a population?",
    options: ["Attack rate", "Incidence rate", "Prevalence", "Mortality rate"],
    answer: 0,
    explanation: "Attack rate measures the proportion of people who become ill after exposure, typically used in outbreak settings."
  },
  {
    question: "Which strategy is used to achieve herd immunity in a population?",
    options: ["Vaccination", "Vector control", "Water chlorination", "Quarantine"],
    answer: 0,
    explanation: "Vaccination induces immunity in a significant portion of the population, reducing disease spread. Vector control, chlorination, and quarantine are other preventive measures."
  },
  {
    question: "Which epidemiological study type measures prevalence at a single point in time to identify associations?",
    options: ["Cross-sectional study", "Case-control study", "Cohort study", "Randomized trial"],
    answer: 0,
    explanation: "Cross-sectional studies measure disease and exposure simultaneously to find correlations. Case-control looks retrospectively, cohort prospectively, randomized trials test interventions."
  },
  {
    question: "Which global event qualifies as a pandemic rather than an epidemic?",
    options: ["Seasonal influenza in one city", "COVID-19 outbreak in multiple continents", "Measles outbreak in a school", "Cholera in a village"],
    answer: 1,
    explanation: "A pandemic spreads across countries/continents (e.g., COVID-19). Localized outbreaks are epidemics or sporadic occurrences."
  },
  {
    question: "Which public health measure is critical to prevent nosocomial infections in hospitals?",
    options: ["Hand hygiene", "Vaccination", "Vector control", "Boiling water"],
    answer: 0,
    explanation: "Hand hygiene prevents transmission of pathogens in hospitals, reducing nosocomial infections. Vaccination prevents community-acquired diseases, vector control targets vectors, boiling water prevents waterborne infections."
  },
  {
    question: "Which method involves using selected sites to monitor disease trends in a population?",
    options: ["Sentinel surveillance", "Active surveillance", "Passive surveillance", "Syndromic surveillance"],
    answer: 0,
    explanation: "Sentinel surveillance uses selected sites or facilities to monitor disease patterns, providing early warning. Active is proactive, passive relies on routine reporting, syndromic tracks symptoms."
  },
  {
    question: "Which environmental factor can increase vector-borne disease transmission?",
    options: ["Climate change", "Vaccination", "Hand hygiene", "Sanitation improvement"],
    answer: 0,
    explanation: "Climate change can expand vector habitats and increase disease transmission. Vaccination, hygiene, and sanitation reduce transmission."
  },
  {
    question: "Which laboratory test is essential for confirming meningococcal outbreaks?",
    options: ["CSF culture for Neisseria meningitidis", "Stool culture for E. coli", "Blood smear for malaria", "KOH mount for fungi"],
    answer: 0,
    explanation: "CSF culture detects Neisseria meningitidis in suspected meningococcal meningitis outbreaks. Other tests target different pathogens."
  }
];

let quizData11 = [
  {
    question: "Which enzyme is commonly used in molecular cloning to cut DNA at specific sequences?",
    options: ["Ligase", "Polymerase", "Restriction endonuclease", "Reverse transcriptase"],
    answer: 2,
    explanation: "Restriction endonucleases recognize specific DNA sequences and cleave them, allowing insertion of genes into vectors. Ligase joins DNA fragments, polymerase synthesizes DNA, reverse transcriptase makes cDNA from RNA."
  },
  {
    question: "Which vector is commonly used to introduce recombinant DNA into bacteria?",
    options: ["Plasmid", "Bacteriophage", "Cosmid", "All of the above"],
    answer: 3,
    explanation: "Plasmids, bacteriophages, and cosmids are all used as vectors for cloning DNA into bacteria depending on size of DNA and application."
  },
  {
    question: "Which microbial host is most widely used for industrial-scale production of recombinant insulin?",
    options: ["Escherichia coli", "Saccharomyces cerevisiae", "Bacillus subtilis", "Aspergillus niger"],
    answer: 0,
    explanation: "Escherichia coli is the primary host for recombinant human insulin production due to rapid growth, well-characterized genetics, and high expression levels."
  },
  {
    question: "Which CRISPR-Cas system component acts as a guide to target specific DNA sequences?",
    options: ["Cas9 nuclease", "sgRNA", "DNA ligase", "Polymerase"],
    answer: 1,
    explanation: "Single guide RNA (sgRNA) directs Cas9 to the complementary DNA sequence for precise genome editing. Cas9 is the nuclease, ligase joins DNA, polymerase synthesizes DNA."
  },
  {
    question: "Which technique allows for amplification of a specific DNA fragment exponentially?",
    options: ["Southern blot", "Northern blot", "PCR", "ELISA"],
    answer: 2,
    explanation: "Polymerase chain reaction (PCR) amplifies specific DNA sequences exponentially. Southern blot detects DNA, Northern blot detects RNA, ELISA detects proteins or antigens."
  },
  {
    question: "Which selectable marker is commonly used in plasmid vectors to identify transformed bacteria?",
    options: ["Antibiotic resistance gene", "Origin of replication", "Promoter sequence", "Terminator sequence"],
    answer: 0,
    explanation: "Antibiotic resistance genes allow selection of transformed bacteria, as only cells with the plasmid can survive on antibiotic-containing media."
  },
  {
    question: "Which technique allows direct sequencing of all microbial DNA in an environmental sample?",
    options: ["Metagenomics", "Sanger sequencing", "PCR", "RFLP"],
    answer: 0,
    explanation: "Metagenomics sequences total DNA from environmental samples to study unculturable microbial communities. Sanger sequencing sequences single genes, PCR amplifies DNA, RFLP detects fragment length polymorphisms."
  },
  {
    question: "Which microbial product is produced using genetically engineered Saccharomyces cerevisiae for bioethanol?",
    options: ["Ethanol", "Citric acid", "Penicillin", "Vitamin B12"],
    answer: 0,
    explanation: "Engineered Saccharomyces cerevisiae ferments sugars into ethanol for biofuel production. Citric acid is from Aspergillus, penicillin from Penicillium, vitamin B12 from Propionibacterium."
  },
  {
    question: "Which gene editing approach allows precise insertion, deletion, or replacement of DNA in microbial genomes?",
    options: ["CRISPR-Cas9", "RNAi", "Transposon mutagenesis", "Random mutagenesis"],
    answer: 0,
    explanation: "CRISPR-Cas9 enables precise genome editing by creating double-strand breaks at targeted locations. RNAi silences gene expression, transposons insert randomly, random mutagenesis is uncontrolled."
  },
  {
    question: "Which microbial biotechnology application produces monoclonal antibodies?",
    options: ["Hybridoma technology", "Fermentation of Lactobacillus", "Ethanol production", "CRISPR-based editing"],
    answer: 0,
    explanation: "Hybridoma technology fuses B-cells with myeloma cells to produce monoclonal antibodies. Fermentation produces metabolites, ethanol production is biofuel, CRISPR edits genomes."
  },
  {
    question: "Which protein is used as a reporter gene in microbial biotechnology for tracking expression?",
    options: ["GFP", "LacZ", "Beta-lactamase", "All of the above"],
    answer: 3,
    explanation: "GFP, LacZ, and beta-lactamase are all used as reporter genes to monitor gene expression, protein localization, or plasmid presence."
  },
  {
    question: "Which microbial system is commonly used for large-scale vaccine production?",
    options: ["Baculovirus-insect cell system", "Escherichia coli", "Yeast expression system", "All of the above"],
    answer: 3,
    explanation: "Baculovirus-insect cell system, E. coli, and yeast are all used for recombinant vaccine production depending on protein complexity and post-translational modifications."
  },
  {
    question: "Which omics approach studies the complete protein set expressed by a microorganism?",
    options: ["Genomics", "Transcriptomics", "Proteomics", "Metabolomics"],
    answer: 2,
    explanation: "Proteomics analyzes all proteins expressed in a cell or organism. Genomics studies DNA, transcriptomics studies RNA, metabolomics studies metabolites."
  },
  {
    question: "Which microbial process is used to synthesize biodegradable plastics like polyhydroxyalkanoates (PHAs)?",
    options: ["Fermentation by Cupriavidus necator", "Methanogenesis", "Alcoholic fermentation", "Denitrification"],
    answer: 0,
    explanation: "Cupriavidus necator and related bacteria produce PHAs intracellularly during fermentation, providing biodegradable plastics."
  },
  {
    question: "Which microbial biotechnology approach enables construction of artificial genetic circuits?",
    options: ["Synthetic biology", "Metagenomics", "Fermentation", "Random mutagenesis"],
    answer: 0,
    explanation: "Synthetic biology designs and constructs artificial genetic circuits in microbes for specific functions. Metagenomics studies natural communities, fermentation produces metabolites, random mutagenesis alters genes randomly."
  },
  {
    question: "Which regulatory element controls the initiation of transcription in recombinant plasmids?",
    options: ["Promoter", "Origin of replication", "Terminator", "Selectable marker"],
    answer: 0,
    explanation: "Promoters are DNA sequences that initiate transcription. Origins replicate DNA, terminators end transcription, selectable markers enable selection."
  },
  {
    question: "Which viral vector is commonly used for gene delivery in microbial and mammalian systems?",
    options: ["Adenovirus", "Retrovirus", "Lentivirus", "All of the above"],
    answer: 3,
    explanation: "Adenovirus, retrovirus, and lentivirus are viral vectors used for gene delivery in microbes and mammalian cells depending on requirements."
  },
  {
    question: "Which method allows introduction of recombinant DNA into bacteria via temporary membrane permeability?",
    options: ["Electroporation", "Conjugation", "Transformation", "Transduction"],
    answer: 0,
    explanation: "Electroporation uses electrical pulses to create temporary pores in bacterial membranes for DNA uptake. Transformation is DNA uptake in general, conjugation involves plasmid transfer via pilus, transduction uses bacteriophages."
  },
  {
    question: "Which microbial biotechnology application involves engineering microbes to produce biofuels from lignocellulosic biomass?",
    options: ["Metabolic engineering of cellulolytic microbes", "Ethanol fermentation by Saccharomyces only", "Methanogenesis only", "All involve only bacteria"],
    answer: 0,
    explanation: "Metabolic engineering of cellulolytic microbes allows degradation of lignocellulose and production of biofuels. Saccharomyces alone cannot utilize cellulose efficiently."
  },
  {
    question: "Which technology allows selective silencing of microbial gene expression without altering DNA sequence?",
    options: ["RNA interference (RNAi)", "CRISPR-Cas9", "Transposon mutagenesis", "Site-directed mutagenesis"],
    answer: 0,
    explanation: "RNAi silences gene expression post-transcriptionally, CRISPR edits DNA, transposons insert randomly, site-directed mutagenesis changes DNA sequence."
  },
  {
    question: "Which bacterial genus is widely used in metabolic engineering for production of amino acids like lysine and glutamate?",
    options: ["Corynebacterium", "Escherichia", "Bacillus", "Pseudomonas"],
    answer: 0,
    explanation: "Corynebacterium glutamicum is industrially engineered to produce large amounts of lysine and glutamate. E. coli is used for recombinant proteins, Bacillus for enzymes, Pseudomonas for biodegradation."
  },
  {
    question: "Which microbial product is obtained via recombinant DNA technology for treatment of diabetes?",
    options: ["Insulin", "Erythropoietin", "Growth hormone", "Interferon"],
    answer: 0,
    explanation: "Recombinant human insulin produced in E. coli or yeast is used to treat diabetes. Erythropoietin treats anemia, growth hormone treats growth deficiency, interferon treats viral infections."
  },
  {
    question: "Which method allows detection of protein-protein interactions in microbial biotechnology?",
    options: ["Yeast two-hybrid system", "PCR", "Southern blot", "Northern blot"],
    answer: 0,
    explanation: "Yeast two-hybrid system detects protein-protein interactions. PCR amplifies DNA, Southern blot detects DNA, Northern blot detects RNA."
  },
  {
    question: "Which microbial biotechnology application uses engineered microbes to degrade environmental pollutants efficiently?",
    options: ["Bioaugmentation", "Fermentation", "Methanogenesis", "Antibiotic production"],
    answer: 0,
    explanation: "Bioaugmentation uses engineered or selected microbes to accelerate degradation of environmental pollutants. Fermentation produces metabolites, methanogenesis produces methane, antibiotic production synthesizes drugs."
  }
];


let quizData12 = [
  {
    question: "Which principle ensures that every batch of microbial product meets predefined quality standards?",
    options: ["Quality Assurance", "Good Laboratory Practice", "Good Manufacturing Practice", "Quality Control"],
    answer: 0,
    explanation: "Quality Assurance (QA) ensures that all processes, procedures, and products meet defined quality standards consistently."
  },
  {
    question: "Which QA component involves routine testing for sterility, purity, and potency of microbial products?",
    options: ["Quality Control", "GLP", "GMP", "SOPs"],
    answer: 0,
    explanation: "Quality Control (QC) involves regular testing and monitoring to ensure product sterility, purity, and potency, while QA encompasses the entire quality system."
  },
  {
    question: "Which regulatory guideline focuses on hazard analysis and critical control points in food microbiology?",
    options: ["ISO 9001", "HACCP", "FDA GMP", "ISO 17025"],
    answer: 1,
    explanation: "HACCP (Hazard Analysis and Critical Control Points) identifies critical points in food production to prevent contamination."
  },
  {
    question: "Which microbial indicator is commonly used to monitor water or food contamination?",
    options: ["E. coli", "Saccharomyces cerevisiae", "Bacillus subtilis", "Trichoderma reesei"],
    answer: 0,
    explanation: "E. coli is widely used as an indicator of fecal contamination and hygiene standards in water and food products."
  },
  {
    question: "Which sterilization method uses high-pressure saturated steam to eliminate all microorganisms including spores?",
    options: ["Autoclaving", "Filtration", "UV irradiation", "Pasteurization"],
    answer: 0,
    explanation: "Autoclaving uses high-pressure steam at 121°C to sterilize equipment and media, killing all microorganisms including spores."
  },
  {
    question: "Which QA standard is specifically used for laboratory competence and testing accuracy?",
    options: ["ISO 17025", "ISO 22000", "ISO 9001", "ISO 13485"],
    answer: 0,
    explanation: "ISO 17025 specifies general requirements for the competence of testing and calibration laboratories."
  },
  {
    question: "Which rapid microbiological method allows detection of pathogens without culture?",
    options: ["PCR", "Viable plate count", "Serial dilution", "Gram staining"],
    answer: 0,
    explanation: "PCR can detect specific DNA sequences of pathogens quickly without the need for culture-based methods."
  },
  {
    question: "Which QA practice ensures reproducibility and reliability of laboratory results?",
    options: ["Good Laboratory Practice (GLP)", "Good Manufacturing Practice (GMP)", "HACCP", "Bioaugmentation"],
    answer: 0,
    explanation: "GLP ensures integrity, reproducibility, and reliability of laboratory studies and experiments."
  },
  {
    question: "Which QA procedure verifies that sterilization processes are effective?",
    options: ["Biological indicators", "Environmental monitoring", "pH measurement", "Temperature log only"],
    answer: 0,
    explanation: "Biological indicators, often spores of Bacillus species, are used to verify that sterilization processes effectively kill all microorganisms."
  },
  {
    question: "Which documentation is essential for traceability in microbiological production?",
    options: ["Batch records", "SOPs", "Environmental logs", "All of the above"],
    answer: 3,
    explanation: "Batch records, SOPs, and environmental monitoring logs together ensure full traceability and compliance in microbial production."
  },
  {
    question: "Which QA approach identifies, evaluates, and controls risks in microbial processes?",
    options: ["Risk-based quality assurance", "Random sampling", "Routine observation", "Metagenomics"],
    answer: 0,
    explanation: "Risk-based QA systematically assesses potential risks and implements control measures to prevent failures."
  },
  {
    question: "Which pathogen is critical to monitor in sterile pharmaceutical products?",
    options: ["Staphylococcus aureus", "Saccharomyces cerevisiae", "Bacillus subtilis", "Trichoderma spp."],
    answer: 0,
    explanation: "Staphylococcus aureus is a significant contaminant in pharmaceuticals and must be monitored to ensure sterility and safety."
  },
  {
    question: "Which microbial QA practice evaluates air and surface contamination in cleanrooms?",
    options: ["Environmental monitoring", "PCR", "Fermentation control", "ELISA"],
    answer: 0,
    explanation: "Environmental monitoring assesses microbial contamination in air, surfaces, and equipment in controlled production environments."
  },
  {
    question: "Which QA document contains step-by-step instructions for laboratory or production procedures?",
    options: ["Standard Operating Procedure (SOP)", "Batch record", "Audit report", "Quality manual"],
    answer: 0,
    explanation: "SOPs provide detailed instructions to ensure consistency and reproducibility in laboratory and industrial procedures."
  },
  {
    question: "Which QA principle ensures that personnel are trained and competent in microbiological procedures?",
    options: ["Training and competency assessment", "Batch record documentation", "Environmental monitoring", "Metabolic engineering"],
    answer: 0,
    explanation: "Training and competency assessment ensures personnel can perform procedures correctly and safely."
  },
  {
    question: "Which QA practice is used to verify the absence of endotoxins in injectable microbial products?",
    options: ["Limulus Amebocyte Lysate (LAL) test", "Gram staining", "PCR", "Fermentation monitoring"],
    answer: 0,
    explanation: "The LAL test detects endotoxins from Gram-negative bacteria to ensure safety of injectable products."
  },
  {
    question: "Which regulatory body oversees pharmaceutical microbial QA in the United States?",
    options: ["FDA", "WHO", "ISO", "CDC"],
    answer: 0,
    explanation: "The FDA regulates quality standards, inspections, and compliance of pharmaceutical microbiology products in the USA."
  },
  {
    question: "Which QA approach involves routine internal and external audits of microbiological laboratories?",
    options: ["Auditing and inspection", "Environmental monitoring", "PCR analysis", "Sterility testing"],
    answer: 0,
    explanation: "Auditing and inspections evaluate compliance with SOPs, regulatory guidelines, and QA systems."
  },
  {
    question: "Which QA concept ensures that all corrective actions for non-conformities are documented and implemented?",
    options: ["Corrective and Preventive Actions (CAPA)", "Environmental monitoring", "Batch record keeping", "PCR testing"],
    answer: 0,
    explanation: "CAPA ensures that deviations are addressed, documented, and prevented from recurring."
  },
  {
    question: "Which method is used to detect yeast and mold contamination in food and pharmaceutical products?",
    options: ["Plating on selective media", "PCR for bacterial DNA", "Gram staining", "Metagenomics only"],
    answer: 0,
    explanation: "Plating on selective media allows enumeration and detection of yeast and molds in products."
  },
  {
    question: "Which microbial QA standard focuses on medical device quality management?",
    options: ["ISO 13485", "ISO 9001", "ISO 22000", "ISO 17025"],
    answer: 0,
    explanation: "ISO 13485 specifies quality management systems requirements for medical devices, including microbiological safety."
  },
  {
    question: "Which QA parameter is essential for ensuring reproducible fermentation processes?",
    options: ["pH, temperature, dissolved oxygen, and microbial inoculum quality", "Airflow only", "Media color", "Time of day"],
    answer: 0,
    explanation: "Monitoring and controlling pH, temperature, dissolved oxygen, and inoculum quality ensures consistent microbial growth and product yield."
  },
  {
    question: "Which QA approach is used to verify that autoclaving effectively sterilized media?",
    options: ["Biological indicator testing", "Gram staining", "Environmental monitoring", "PCR detection"],
    answer: 0,
    explanation: "Biological indicators, usually resistant spores, verify the effectiveness of autoclaving."
  },
  {
    question: "Which QA practice ensures consistency in microbial assay results across multiple laboratories?",
    options: ["Inter-laboratory validation", "Single lab testing", "Random sampling", "Pure culture isolation only"],
    answer: 0,
    explanation: "Inter-laboratory validation ensures reproducibility, accuracy, and comparability of results among different labs."
  },
  {
    question: "Which rapid detection method in QA can quantify viable microorganisms without culturing?",
    options: ["Flow cytometry", "Gram staining", "Serial dilution", "Metagenomics"],
    answer: 0,
    explanation: "Flow cytometry allows rapid quantification and analysis of live microbial cells without the need for traditional culture."
  }
];


/*************************************************
 * 2️⃣ AUTO-COLLECT ALL QUIZ DATA
 *************************************************/

const allQuizData = [
  ...quizData1,
  ...quizData2,
  ...quizData3,
  ...quizData4,
  ...quizData5,
  ...quizData6,
  ...quizData7,
  ...quizData8,
  ...quizData9,
  ...quizData10,
  ...quizData11,
  ...quizData12
];