#!/usr/bin/env python3
"""
Enhanced ScType Database Creator
Compiles markers from CellMarker 2.0, PanglaoDB, Human Cell Atlas, and literature
Creates comprehensive database for both human and mouse cells
"""

import pandas as pd
from pathlib import Path

def create_enhanced_database():
    """Create comprehensive marker database from multiple sources"""

    # Initialize list to collect all cell type data
    all_cells = []

    #============================================================================#
    # IMMUNE SYSTEM - Extensively updated from CellMarker 2.0 & literature
    #============================================================================#

    immune_cells = [
        # T cells
        {"tissueType": "Immune system", "cellName": "T cells", "shortName": "T",
         "geneSymbolmore1": "CD3D,CD3E,CD3G,CD2,CD5,CD7,TRAC,TRBC1,TRBC2",
         "geneSymbolmore2": "CD19,CD20,MS4A1,CD79A,CD14,FCGR3A"},

        {"tissueType": "Immune system", "cellName": "CD4+ T cells", "shortName": "CD4T",
         "geneSymbolmore1": "CD4,CD3D,CD3E,IL7R,MAL,LTB",
         "geneSymbolmore2": "CD8A,CD8B,NKG7,GNLY,NCAM1,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "CD8+ T cells", "shortName": "CD8T",
         "geneSymbolmore1": "CD8A,CD8B,CD3D,CD3E,GZMA,GZMB,PRF1,NKG7",
         "geneSymbolmore2": "CD4,IL7R,CD19,CD14,FCGR3A"},

        {"tissueType": "Immune system", "cellName": "Naive CD4+ T cells", "shortName": "NaiveCD4T",
         "geneSymbolmore1": "CCR7,SELL,LEF1,TCF7,CD4,IL7R",
         "geneSymbolmore2": "CD44,CD69,CD8A,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "Naive CD8+ T cells", "shortName": "NaiveCD8T",
         "geneSymbolmore1": "CCR7,SELL,LEF1,TCF7,CD8A,IL7R",
         "geneSymbolmore2": "CD44,CD69,CD4,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "Memory CD4+ T cells", "shortName": "MemCD4T",
         "geneSymbolmore1": "CD4,IL7R,CD44,CD69,LTB,AQP3",
         "geneSymbolmore2": "CCR7,SELL,CD8A,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "Memory CD8+ T cells", "shortName": "MemCD8T",
         "geneSymbolmore1": "CD8A,CD8B,CD44,GZMK,CCL5,DUSP2",
         "geneSymbolmore2": "CCR7,SELL,CD4,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "Regulatory T cells", "shortName": "Treg",
         "geneSymbolmore1": "FOXP3,IL2RA,CTLA4,IKZF2,CD4,CD25,TNFRSF18",
         "geneSymbolmore2": "CD8A,IL17A,IFNG,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "Th1 cells", "shortName": "Th1",
         "geneSymbolmore1": "IFNG,TBX21,CXCR3,CCR5,CD4,IL12RB2",
         "geneSymbolmore2": "IL4,GATA3,IL17A,RORC,CD8A,CD19"},

        {"tissueType": "Immune system", "cellName": "Th2 cells", "shortName": "Th2",
         "geneSymbolmore1": "GATA3,IL4,IL5,IL13,CCR4,CD4",
         "geneSymbolmore2": "IFNG,TBX21,IL17A,RORC,CD8A,CD19"},

        {"tissueType": "Immune system", "cellName": "Th17 cells", "shortName": "Th17",
         "geneSymbolmore1": "IL17A,IL17F,RORC,CCR6,CD4,IL23R",
         "geneSymbolmore2": "IFNG,TBX21,GATA3,IL4,CD8A,CD19"},

        {"tissueType": "Immune system", "cellName": "Follicular helper T cells", "shortName": "Tfh",
         "geneSymbolmore1": "CXCR5,BCL6,PDCD1,ICOS,CD4,IL21",
         "geneSymbolmore2": "CD8A,FOXP3,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "Gamma delta T cells", "shortName": "gdT",
         "geneSymbolmore1": "TRDC,TRDV1,TRGC1,TRGC2,TRGV9",
         "geneSymbolmore2": "CD4,CD8A,TRAC,TRBC1,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "MAIT cells", "shortName": "MAIT",
         "geneSymbolmore1": "SLC4A10,KLRB1,ZBTB16,TRAV1-2,NCR3",
         "geneSymbolmore2": "CD4,CD19,CD14"},

        # B cells
        {"tissueType": "Immune system", "cellName": "B cells", "shortName": "B",
         "geneSymbolmore1": "CD19,MS4A1,CD79A,CD79B,CD22",
         "geneSymbolmore2": "CD3D,CD3E,CD14,FCGR3A,NKG7"},

        {"tissueType": "Immune system", "cellName": "Naive B cells", "shortName": "NaiveB",
         "geneSymbolmore1": "CD19,MS4A1,TCL1A,IGHD,IGHM,IL4R",
         "geneSymbolmore2": "CD27,CD38,MZB1,XBP1,CD3D,CD14"},

        {"tissueType": "Immune system", "cellName": "Memory B cells", "shortName": "MemB",
         "geneSymbolmore1": "CD19,MS4A1,CD27,TNFRSF13B,AIM2",
         "geneSymbolmore2": "IGHD,TCL1A,CD38,MZB1,CD3D,CD14"},

        {"tissueType": "Immune system", "cellName": "Plasma cells", "shortName": "Plasma",
         "geneSymbolmore1": "MZB1,XBP1,SDC1,CD27,CD38,JCHAIN,IGHG1,IGHA1",
         "geneSymbolmore2": "CD19,MS4A1,CD79A,CD3D,CD14"},

        {"tissueType": "Immune system", "cellName": "Plasmablasts", "shortName": "Plasmablast",
         "geneSymbolmore1": "XBP1,PRDM1,CD27,CD38,MZB1,JCHAIN",
         "geneSymbolmore2": "MS4A1,TCL1A,CD3D,CD14"},

        # NK cells
        {"tissueType": "Immune system", "cellName": "NK cells", "shortName": "NK",
         "geneSymbolmore1": "NKG7,GNLY,KLRD1,KLRF1,NCR1,NCAM1,FCGR3A",
         "geneSymbolmore2": "CD3D,CD3E,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "CD56bright NK cells", "shortName": "CD56brightNK",
         "geneSymbolmore1": "NCAM1,XCL1,XCL2,SELL,CD7,KLRC1",
         "geneSymbolmore2": "FCGR3A,GZMB,PRF1,CD3D,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "CD56dim NK cells", "shortName": "CD56dimNK",
         "geneSymbolmore1": "FCGR3A,GZMB,PRF1,KLRF1,KLRD1,CX3CR1",
         "geneSymbolmore2": "XCL1,XCL2,SELL,CD3D,CD19,CD14"},

        # Monocytes and Macrophages
        {"tissueType": "Immune system", "cellName": "Monocytes", "shortName": "Mono",
         "geneSymbolmore1": "CD14,LYZ,S100A8,S100A9,VCAN,FCN1",
         "geneSymbolmore2": "CD3D,CD19,NCAM1,NKG7"},

        {"tissueType": "Immune system", "cellName": "Classical monocytes", "shortName": "cMono",
         "geneSymbolmore1": "CD14,S100A8,S100A9,LYZ,VCAN,FCN1",
         "geneSymbolmore2": "FCGR3A,MS4A7,CD16,CD3D,CD19"},

        {"tissueType": "Immune system", "cellName": "Non-classical monocytes", "shortName": "ncMono",
         "geneSymbolmore1": "FCGR3A,MS4A7,CDKN1C,CX3CR1,PTPRC",
         "geneSymbolmore2": "CD14,S100A8,S100A9,CD3D,CD19"},

        {"tissueType": "Immune system", "cellName": "Intermediate monocytes", "shortName": "intMono",
         "geneSymbolmore1": "CD14,FCGR3A,HLA-DRA,CD74,CST3",
         "geneSymbolmore2": "CD3D,CD19,NKG7"},

        {"tissueType": "Immune system", "cellName": "Macrophages", "shortName": "Macro",
         "geneSymbolmore1": "CD68,CD163,MSR1,MRC1,C1QA,C1QB,C1QC",
         "geneSymbolmore2": "CD3D,CD19,NKG7,CD79A"},

        {"tissueType": "Immune system", "cellName": "M1 Macrophages", "shortName": "M1Macro",
         "geneSymbolmore1": "CD68,CD80,CD86,IL1B,TNF,CXCL10,IL12A",
         "geneSymbolmore2": "CD163,MRC1,IL10,CD3D,CD19"},

        {"tissueType": "Immune system", "cellName": "M2 Macrophages", "shortName": "M2Macro",
         "geneSymbolmore1": "CD68,CD163,MRC1,MSR1,IL10,TGFB1,ARG1",
         "geneSymbolmore2": "IL1B,TNF,CXCL10,CD3D,CD19"},

        # Dendritic cells
        {"tissueType": "Immune system", "cellName": "Dendritic cells", "shortName": "DC",
         "geneSymbolmore1": "FCER1A,CST3,CD1C,CLEC10A,HLA-DRA,HLA-DRB1",
         "geneSymbolmore2": "CD3D,CD19,CD14,NKG7"},

        {"tissueType": "Immune system", "cellName": "Conventional DC1", "shortName": "cDC1",
         "geneSymbolmore1": "CLEC9A,XCR1,BATF3,IRF8,ID2,CADM1",
         "geneSymbolmore2": "CD1C,FCER1A,CD14,CD3D,CD19"},

        {"tissueType": "Immune system", "cellName": "Conventional DC2", "shortName": "cDC2",
         "geneSymbolmore1": "CD1C,FCER1A,CLEC10A,CD1E,IRF4",
         "geneSymbolmore2": "CLEC9A,XCR1,CD141,CD14,CD3D,CD19"},

        {"tissueType": "Immune system", "cellName": "Plasmacytoid DC", "shortName": "pDC",
         "geneSymbolmore1": "IL3RA,CLEC4C,LILRA4,IRF7,IRF8,TCF4",
         "geneSymbolmore2": "CD1C,CLEC9A,CD14,CD3D,CD19"},

        {"tissueType": "Immune system", "cellName": "Activated DC", "shortName": "actDC",
         "geneSymbolmore1": "CCR7,LAMP3,CCL19,CCL22,FSCN1,CD86",
         "geneSymbolmore2": "CD14,CD3D,CD19"},

        # Granulocytes
        {"tissueType": "Immune system", "cellName": "Neutrophils", "shortName": "Neutro",
         "geneSymbolmore1": "FCGR3B,CSF3R,FPR1,S100A12,CXCR2,CD177",
         "geneSymbolmore2": "CD3D,CD19,CD14,NCAM1"},

        {"tissueType": "Immune system", "cellName": "Eosinophils", "shortName": "Eos",
         "geneSymbolmore1": "CCR3,SIGLEC8,IL5RA,CLC,PRG2,EPX",
         "geneSymbolmore2": "CD3D,CD19,CD14,NCAM1"},

        {"tissueType": "Immune system", "cellName": "Basophils", "shortName": "Baso",
         "geneSymbolmore1": "MS4A2,CPA3,TPSAB1,TPSB2,HDC,IL3RA",
         "geneSymbolmore2": "CD3D,CD19,CD14,NCAM1"},

        {"tissueType": "Immune system", "cellName": "Mast cells", "shortName": "Mast",
         "geneSymbolmore1": "TPSAB1,TPSB2,CPA3,KIT,MS4A2,HDC",
         "geneSymbolmore2": "CD3D,CD19,CD14,NCAM1,PTPRC"},

        # ILCs
        {"tissueType": "Immune system", "cellName": "ILC1", "shortName": "ILC1",
         "geneSymbolmore1": "TBX21,IFNG,NCR1,IL12RB2,EOMES",
         "geneSymbolmore2": "GATA3,RORC,CD3D,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "ILC2", "shortName": "ILC2",
         "geneSymbolmore1": "GATA3,IL5,IL13,PTGDR2,IL1RL1",
         "geneSymbolmore2": "TBX21,RORC,CD3D,CD19,CD14"},

        {"tissueType": "Immune system", "cellName": "ILC3", "shortName": "ILC3",
         "geneSymbolmore1": "RORC,IL22,IL23R,KIT,CCR6",
         "geneSymbolmore2": "TBX21,GATA3,CD3D,CD19,CD14"},

        # HSCs and blood cells
        {"tissueType": "Immune system", "cellName": "Hematopoietic stem cells", "shortName": "HSC",
         "geneSymbolmore1": "CD34,THY1,KIT,CD133,PROM1,ENG",
         "geneSymbolmore2": "CD3D,CD19,CD14,CD38,HLA-DR"},

        {"tissueType": "Immune system", "cellName": "Erythrocytes", "shortName": "RBC",
         "geneSymbolmore1": "HBA1,HBA2,HBB,HBD,GYPA,SLC4A1",
         "geneSymbolmore2": "CD3D,CD19,CD14,PTPRC"},

        {"tissueType": "Immune system", "cellName": "Platelets", "shortName": "Platelet",
         "geneSymbolmore1": "PPBP,PF4,GP1BA,GP9,ITGA2B,TUBB1",
         "geneSymbolmore2": "CD3D,CD19,CD14,PTPRC"},
    ]

    all_cells.extend(immune_cells)

    #============================================================================#
    # BRAIN - Updated from literature and Human Brain Atlas
    #============================================================================#

    brain_cells = [
        {"tissueType": "Brain", "cellName": "Neurons", "shortName": "Neuron",
         "geneSymbolmore1": "RBFOX3,SYT1,SNAP25,SLC17A7,SLC17A6,STMN2,MAP2",
         "geneSymbolmore2": "GFAP,AQP4,MBP,MOG,CD68,PTPRC"},

        {"tissueType": "Brain", "cellName": "Excitatory neurons", "shortName": "ExNeuron",
         "geneSymbolmore1": "SLC17A7,SLC17A6,SATB2,TBR1,NEUROD6,CUX2",
         "geneSymbolmore2": "GAD1,GAD2,SLC32A1,GFAP,MBP"},

        {"tissueType": "Brain", "cellName": "Inhibitory neurons", "shortName": "InNeuron",
         "geneSymbolmore1": "GAD1,GAD2,SLC32A1,DLX1,DLX2,LHX6",
         "geneSymbolmore2": "SLC17A7,SLC17A6,GFAP,MBP"},

        {"tissueType": "Brain", "cellName": "Dopaminergic neurons", "shortName": "DANeuron",
         "geneSymbolmore1": "TH,DDC,SLC6A3,DRD2,EN1,LMX1B",
         "geneSymbolmore2": "SLC17A7,GAD1,GFAP,MBP"},

        {"tissueType": "Brain", "cellName": "GABAergic neurons", "shortName": "GABANeuron",
         "geneSymbolmore1": "GAD1,GAD2,SLC32A1,DLX1,SST,PVALB,VIP",
         "geneSymbolmore2": "SLC17A7,SLC17A6,GFAP,MBP"},

        {"tissueType": "Brain", "cellName": "Serotonergic neurons", "shortName": "5HTNeuron",
         "geneSymbolmore1": "TPH2,SLC6A4,DDC,FEV,GATA3",
         "geneSymbolmore2": "TH,DBH,GFAP,MBP"},

        {"tissueType": "Brain", "cellName": "Cholinergic neurons", "shortName": "ChAT",
         "geneSymbolmore1": "CHAT,SLC5A7,SLC18A3,ISL1",
         "geneSymbolmore2": "TH,DBH,GFAP,MBP"},

        {"tissueType": "Brain", "cellName": "Astrocytes", "shortName": "Astro",
         "geneSymbolmore1": "GFAP,AQP4,SLC1A3,SLC1A2,ALDH1L1,S100B",
         "geneSymbolmore2": "RBFOX3,MBP,MOG,OLIG2,CD68"},

        {"tissueType": "Brain", "cellName": "Protoplasmic astrocytes", "shortName": "ProtoAstro",
         "geneSymbolmore1": "GFAP,AQP4,SLC1A2,GJA1,ID3,MFGE8",
         "geneSymbolmore2": "RBFOX3,MBP,OLIG2"},

        {"tissueType": "Brain", "cellName": "Fibrous astrocytes", "shortName": "FibroAstro",
         "geneSymbolmore1": "GFAP,AQP4,SLC1A3,CD44,VIM,CRYAB",
         "geneSymbolmore2": "RBFOX3,MBP,OLIG2"},

        {"tissueType": "Brain", "cellName": "Oligodendrocytes", "shortName": "Oligo",
         "geneSymbolmore1": "MBP,MOG,MAG,PLP1,MOBP,OPALIN",
         "geneSymbolmore2": "GFAP,RBFOX3,CD68,PDGFRA"},

        {"tissueType": "Brain", "cellName": "Oligodendrocyte precursor cells", "shortName": "OPC",
         "geneSymbolmore1": "PDGFRA,CSPG4,VCAN,OLIG1,OLIG2,SOX10",
         "geneSymbolmore2": "MBP,MOG,PLP1,GFAP,RBFOX3"},

        {"tissueType": "Brain", "cellName": "Microglia", "shortName": "Micro",
         "geneSymbolmore1": "TMEM119,P2RY12,CX3CR1,AIF1,CSF1R,PTPRC",
         "geneSymbolmore2": "GFAP,MBP,RBFOX3,PDGFRA"},

        {"tissueType": "Brain", "cellName": "Ependymal cells", "shortName": "Ependymal",
         "geneSymbolmore1": "FOXJ1,CCDC153,DNAH5,RSPH1,PIFO",
         "geneSymbolmore2": "GFAP,MBP,RBFOX3,AQP4"},

        {"tissueType": "Brain", "cellName": "Endothelial cells", "shortName": "Endo",
         "geneSymbolmore1": "PECAM1,VWF,CDH5,FLT1,CLDN5,ENG",
         "geneSymbolmore2": "GFAP,MBP,RBFOX3,AQP4"},

        {"tissueType": "Brain", "cellName": "Pericytes", "shortName": "Pericyte",
         "geneSymbolmore1": "PDGFRB,ANPEP,RGS5,ACTA2,CSPG4",
         "geneSymbolmore2": "PECAM1,GFAP,MBP,RBFOX3"},

        {"tissueType": "Brain", "cellName": "Vascular smooth muscle cells", "shortName": "VSMC",
         "geneSymbolmore1": "ACTA2,MYH11,TAGLN,CNN1,PDGFRB",
         "geneSymbolmore2": "PECAM1,GFAP,MBP,RBFOX3"},
    ]

    all_cells.extend(brain_cells)

    #============================================================================#
    # LIVER - From Human Liver Cell Atlas
    #============================================================================#

    liver_cells = [
        {"tissueType": "Liver", "cellName": "Hepatocytes", "shortName": "Hepato",
         "geneSymbolmore1": "ALB,APOA1,APOB,TTR,APOE,HNF4A,TF,ASGR1",
         "geneSymbolmore2": "KRT7,KRT19,EPCAM,CD68,PECAM1,PTPRC"},

        {"tissueType": "Liver", "cellName": "Periportal hepatocytes", "shortName": "PP_Hepato",
         "geneSymbolmore1": "ALB,CYP2E1,CYP2F1,PCK1,HAL,ASS1",
         "geneSymbolmore2": "CYP2E1,GLUL,KRT19,EPCAM"},

        {"tissueType": "Liver", "cellName": "Pericentral hepatocytes", "shortName": "PC_Hepato",
         "geneSymbolmore1": "ALB,CYP2E1,GLUL,OAT,CYP1A2,CYP2A6",
         "geneSymbolmore2": "PCK1,ASS1,KRT19,EPCAM"},

        {"tissueType": "Liver", "cellName": "Cholangiocytes", "shortName": "Chol",
         "geneSymbolmore1": "KRT7,KRT19,EPCAM,CFTR,SOX9,CLDN4",
         "geneSymbolmore2": "ALB,HNF4A,CD68,PECAM1"},

        {"tissueType": "Liver", "cellName": "Hepatic stellate cells", "shortName": "HSC",
         "geneSymbolmore1": "ACTA2,COL1A1,COL1A2,VIM,PDGFRB,LRAT,HAND2",
         "geneSymbolmore2": "ALB,KRT19,PECAM1,CD68,EPCAM"},

        {"tissueType": "Liver", "cellName": "Kupffer cells", "shortName": "Kupffer",
         "geneSymbolmore1": "CD68,CD163,MARCO,VSIG4,TIMD4,CLEC4F",
         "geneSymbolmore2": "ALB,KRT19,PECAM1,EPCAM,ACTA2"},

        {"tissueType": "Liver", "cellName": "Liver sinusoidal endothelial cells", "shortName": "LSEC",
         "geneSymbolmore1": "PECAM1,VWF,CDH5,FCN2,FCN3,CLEC4G,STAB2",
         "geneSymbolmore2": "ALB,KRT19,CD68,ACTA2,EPCAM"},

        {"tissueType": "Liver", "cellName": "Hepatic progenitor cells", "shortName": "HepProg",
         "geneSymbolmore1": "EPCAM,KRT19,SOX9,PROM1,CD133,ANXA4",
         "geneSymbolmore2": "ALB,HNF4A,CD68,PECAM1"},
    ]

    all_cells.extend(liver_cells)

    #============================================================================#
    # PANCREAS - From literature and pancreatic atlases
    #============================================================================#

    pancreas_cells = [
        {"tissueType": "Pancreas", "cellName": "Beta cells", "shortName": "Beta",
         "geneSymbolmore1": "INS,IAPP,MAFA,NKX6-1,PDX1,SLC30A8",
         "geneSymbolmore2": "GCG,SST,PPY,AMY2A,KRT19"},

        {"tissueType": "Pancreas", "cellName": "Alpha cells", "shortName": "Alpha",
         "geneSymbolmore1": "GCG,ARX,IRX2,TTR,GC,PCSK1",
         "geneSymbolmore2": "INS,SST,PPY,AMY2A,KRT19"},

        {"tissueType": "Pancreas", "cellName": "Delta cells", "shortName": "Delta",
         "geneSymbolmore1": "SST,HHEX,LEPR,PRG4,RBP4",
         "geneSymbolmore2": "INS,GCG,PPY,AMY2A,KRT19"},

        {"tissueType": "Pancreas", "cellName": "PP cells", "shortName": "PP",
         "geneSymbolmore1": "PPY,MEIS2,SERTM1,CARTPT",
         "geneSymbolmore2": "INS,GCG,SST,AMY2A,KRT19"},

        {"tissueType": "Pancreas", "cellName": "Epsilon cells", "shortName": "Epsilon",
         "geneSymbolmore1": "GHRL,NEUROD1,MLN,CHGA",
         "geneSymbolmore2": "INS,GCG,SST,PPY,AMY2A"},

        {"tissueType": "Pancreas", "cellName": "Acinar cells", "shortName": "Acinar",
         "geneSymbolmore1": "AMY2A,PRSS1,CPA1,CPA2,REG1A,CTRB1",
         "geneSymbolmore2": "INS,GCG,SST,KRT19,EPCAM"},

        {"tissueType": "Pancreas", "cellName": "Ductal cells", "shortName": "Ductal",
         "geneSymbolmore1": "KRT19,CFTR,SOX9,MMP7,SPP1,AMBP",
         "geneSymbolmore2": "INS,GCG,SST,AMY2A,PRSS1"},

        {"tissueType": "Pancreas", "cellName": "Pancreatic stellate cells", "shortName": "PSC",
         "geneSymbolmore1": "COL1A1,COL1A2,ACTA2,VIM,PDGFRB,RGS5",
         "geneSymbolmore2": "INS,KRT19,AMY2A,PECAM1"},

        {"tissueType": "Pancreas", "cellName": "Endothelial cells", "shortName": "Endo",
         "geneSymbolmore1": "PECAM1,VWF,CDH5,PLVAP,ENG,FLT1",
         "geneSymbolmore2": "INS,KRT19,AMY2A,ACTA2"},
    ]

    all_cells.extend(pancreas_cells)

    #============================================================================#
    # KIDNEY - From Human Kidney Cell Atlas
    #============================================================================#

    kidney_cells = [
        {"tissueType": "Kidney", "cellName": "Podocytes", "shortName": "Podo",
         "geneSymbolmore1": "NPHS1,NPHS2,PODXL,MAGI2,SYNPO,WT1",
         "geneSymbolmore2": "CUBN,SLC34A1,AQP2,SLC12A1,PECAM1"},

        {"tissueType": "Kidney", "cellName": "Proximal tubule cells", "shortName": "PT",
         "geneSymbolmore1": "LRP2,CUBN,SLC34A1,SLC5A2,SLC22A6,GATM",
         "geneSymbolmore2": "NPHS1,AQP2,SLC12A1,PECAM1"},

        {"tissueType": "Kidney", "cellName": "Loop of Henle", "shortName": "LOH",
         "geneSymbolmore1": "SLC12A1,UMOD,CLDN16,CLDN19,PTGER3",
         "geneSymbolmore2": "CUBN,AQP2,NPHS1,PECAM1"},

        {"tissueType": "Kidney", "cellName": "Distal convoluted tubule", "shortName": "DCT",
         "geneSymbolmore1": "SLC12A3,TRPM6,CALB1,PVALB,PCP4",
         "geneSymbolmore2": "CUBN,AQP2,SLC12A1,NPHS1"},

        {"tissueType": "Kidney", "cellName": "Collecting duct principal cells", "shortName": "PC",
         "geneSymbolmore1": "AQP2,AQP3,SCNN1A,SCNN1B,FXYD4",
         "geneSymbolmore2": "CUBN,NPHS1,SLC12A1,ATP6V1B1"},

        {"tissueType": "Kidney", "cellName": "Collecting duct intercalated cells", "shortName": "IC",
         "geneSymbolmore1": "ATP6V1B1,ATP6V0D2,SLC4A1,SLC26A4,FOXI1",
         "geneSymbolmore2": "AQP2,CUBN,NPHS1,SLC12A1"},

        {"tissueType": "Kidney", "cellName": "Mesangial cells", "shortName": "Mesangial",
         "geneSymbolmore1": "PDGFRB,ITGA8,IGFBP5,REN,COL4A1,COL4A2",
         "geneSymbolmore2": "NPHS1,CUBN,PECAM1,AQP2"},

        {"tissueType": "Kidney", "cellName": "Endothelial cells", "shortName": "Endo",
         "geneSymbolmore1": "PECAM1,VWF,CDH5,FLT1,PLVAP,ENG",
         "geneSymbolmore2": "NPHS1,CUBN,AQP2,ACTA2"},

        {"tissueType": "Kidney", "cellName": "Glomerular endothelial cells", "shortName": "GEC",
         "geneSymbolmore1": "PECAM1,PLVAP,PODXL,EMCN,FLT1,EHHADH",
         "geneSymbolmore2": "NPHS1,CUBN,AQP2,ACTA2"},

        {"tissueType": "Kidney", "cellName": "Parietal epithelial cells", "shortName": "PEC",
         "geneSymbolmore1": "PAX2,PAX8,CLDN1,CDH1,KRT8,KRT18",
         "geneSymbolmore2": "NPHS1,CUBN,AQP2,PECAM1"},
    ]

    all_cells.extend(kidney_cells)

    #============================================================================#
    # LUNG - From Human Lung Cell Atlas
    #============================================================================#

    lung_cells = [
        {"tissueType": "Lung", "cellName": "AT1 cells", "shortName": "AT1",
         "geneSymbolmore1": "AGER,PDPN,CAV1,AQP5,HOPX,RTKN2",
         "geneSymbolmore2": "SFTPC,SFTPB,SCGB1A1,FOXJ1,MUC5AC"},

        {"tissueType": "Lung", "cellName": "AT2 cells", "shortName": "AT2",
         "geneSymbolmore1": "SFTPC,SFTPB,SFTPA1,SFTPA2,SFTPD,ABCA3,LPCAT1",
         "geneSymbolmore2": "AGER,SCGB1A1,FOXJ1,MUC5AC,PDPN"},

        {"tissueType": "Lung", "cellName": "Club cells", "shortName": "Club",
         "geneSymbolmore1": "SCGB1A1,SCGB3A2,CC10,CYP2F2,MUC5B,LYPD2",
         "geneSymbolmore2": "SFTPC,FOXJ1,MUC5AC,AGER,KRT5"},

        {"tissueType": "Lung", "cellName": "Ciliated cells", "shortName": "Ciliated",
         "geneSymbolmore1": "FOXJ1,CCDC153,TUBA1A,DNAH5,RSPH1,PIFO",
         "geneSymbolmore2": "SCGB1A1,SFTPC,MUC5AC,KRT5,AGER"},

        {"tissueType": "Lung", "cellName": "Goblet cells", "shortName": "Goblet",
         "geneSymbolmore1": "MUC5AC,MUC5B,TFF3,SPDEF,AGR2",
         "geneSymbolmore2": "SCGB1A1,SFTPC,FOXJ1,KRT5,AGER"},

        {"tissueType": "Lung", "cellName": "Basal cells", "shortName": "Basal",
         "geneSymbolmore1": "KRT5,KRT14,TP63,NGFR,DLK2,KRT17",
         "geneSymbolmore2": "SCGB1A1,SFTPC,FOXJ1,MUC5AC,AGER"},

        {"tissueType": "Lung", "cellName": "Neuroendocrine cells", "shortName": "NE",
         "geneSymbolmore1": "CHGA,CHGB,SYP,NCAM1,ASCL1,GRP",
         "geneSymbolmore2": "SCGB1A1,SFTPC,KRT5,AGER"},

        {"tissueType": "Lung", "cellName": "Alveolar macrophages", "shortName": "AlvMacro",
         "geneSymbolmore1": "MARCO,MSR1,FABP4,CD68,PPARG,MERTK",
         "geneSymbolmore2": "SFTPC,KRT5,SCGB1A1,CD3D"},

        {"tissueType": "Lung", "cellName": "Endothelial cells", "shortName": "Endo",
         "geneSymbolmore1": "PECAM1,VWF,CDH5,ENG,FLT1,CLDN5",
         "geneSymbolmore2": "SFTPC,KRT5,CD68,ACTA2"},

        {"tissueType": "Lung", "cellName": "Fibroblasts", "shortName": "Fibro",
         "geneSymbolmore1": "COL1A1,COL1A2,COL3A1,DCN,LUM,PDGFRA",
         "geneSymbolmore2": "PECAM1,SFTPC,KRT5,CD68"},

        {"tissueType": "Lung", "cellName": "Smooth muscle cells", "shortName": "SMC",
         "geneSymbolmore1": "ACTA2,MYH11,TAGLN,CNN1,PDGFRB",
         "geneSymbolmore2": "PECAM1,SFTPC,KRT5,CD68"},
    ]

    all_cells.extend(lung_cells)

    #============================================================================#
    # HEART - From Human Heart Cell Atlas
    #============================================================================#

    heart_cells = [
        {"tissueType": "Heart", "cellName": "Cardiomyocytes", "shortName": "CM",
         "geneSymbolmore1": "TNNT2,MYH6,MYH7,TNNI3,ACTC1,MYL2,MYL7",
         "geneSymbolmore2": "PECAM1,CD68,PDGFRA,COL1A1"},

        {"tissueType": "Heart", "cellName": "Atrial cardiomyocytes", "shortName": "ACM",
         "geneSymbolmore1": "NPPA,NPPB,MYL4,MYL7,TNNT2,GJA5",
         "geneSymbolmore2": "MYH7,IRX4,PECAM1,CD68"},

        {"tissueType": "Heart", "cellName": "Ventricular cardiomyocytes", "shortName": "VCM",
         "geneSymbolmore1": "MYH7,IRX4,HEY2,MYL2,TNNT2,GJA1",
         "geneSymbolmore2": "NPPA,MYL7,MYL4,PECAM1,CD68"},

        {"tissueType": "Heart", "cellName": "Cardiac fibroblasts", "shortName": "CF",
         "geneSymbolmore1": "COL1A1,COL1A2,COL3A1,DCN,PDGFRA,THY1",
         "geneSymbolmore2": "TNNT2,PECAM1,CD68,ACTA2"},

        {"tissueType": "Heart", "cellName": "Endothelial cells", "shortName": "Endo",
         "geneSymbolmore1": "PECAM1,VWF,CDH5,KDR,FLT1,ENG",
         "geneSymbolmore2": "TNNT2,CD68,COL1A1,ACTA2"},

        {"tissueType": "Heart", "cellName": "Pericytes", "shortName": "Peri",
         "geneSymbolmore1": "RGS5,PDGFRB,CSPG4,ANPEP,NOTCH3",
         "geneSymbolmore2": "PECAM1,TNNT2,CD68,COL1A1"},

        {"tissueType": "Heart", "cellName": "Smooth muscle cells", "shortName": "SMC",
         "geneSymbolmore1": "ACTA2,MYH11,TAGLN,CNN1,PDGFRB",
         "geneSymbolmore2": "PECAM1,TNNT2,CD68,COL1A1"},

        {"tissueType": "Heart", "cellName": "Adipocytes", "shortName": "Adipo",
         "geneSymbolmore1": "ADIPOQ,LEP,PLIN1,FABP4,LPL,PPARG",
         "geneSymbolmore2": "TNNT2,PECAM1,CD68,COL1A1"},
    ]

    all_cells.extend(heart_cells)

    #============================================================================#
    # INTESTINE
    #============================================================================#

    intestine_cells = [
        {"tissueType": "Intestine", "cellName": "Enterocytes", "shortName": "Entero",
         "geneSymbolmore1": "ALPI,APOA1,FABP2,SI,FABP1,VIL1",
         "geneSymbolmore2": "MUC2,LYZ,CHGA,OLFM4,PECAM1"},

        {"tissueType": "Intestine", "cellName": "Goblet cells", "shortName": "Goblet",
         "geneSymbolmore1": "MUC2,TFF3,FCGBP,ZG16,AGR2",
         "geneSymbolmore2": "ALPI,LYZ,CHGA,KRT20,PECAM1"},

        {"tissueType": "Intestine", "cellName": "Paneth cells", "shortName": "Paneth",
         "geneSymbolmore1": "LYZ,DEFA5,DEFA6,PLA2G2A,PRSS2",
         "geneSymbolmore2": "MUC2,CHGA,OLFM4,PECAM1"},

        {"tissueType": "Intestine", "cellName": "Enteroendocrine cells", "shortName": "EE",
         "geneSymbolmore1": "CHGA,CHGB,TAC1,TPH1,NEUROD1,GCG",
         "geneSymbolmore2": "MUC2,LYZ,ALPI,OLFM4,PECAM1"},

        {"tissueType": "Intestine", "cellName": "Intestinal stem cells", "shortName": "ISC",
         "geneSymbolmore1": "LGR5,OLFM4,ASCL2,SOX9,SMOC2",
         "geneSymbolmore2": "MUC2,LYZ,CHGA,ALPI,PECAM1"},

        {"tissueType": "Intestine", "cellName": "Tuft cells", "shortName": "Tuft",
         "geneSymbolmore1": "DCLK1,TRPM5,GFI1B,IL25,POU2F3",
         "geneSymbolmore2": "MUC2,LYZ,ALPI,PECAM1"},
    ]

    all_cells.extend(intestine_cells)

    #============================================================================#
    # MUSCLE
    #============================================================================#

    muscle_cells = [
        {"tissueType": "Muscle", "cellName": "Skeletal muscle cells", "shortName": "SkMuscle",
         "geneSymbolmore1": "MYH1,MYH2,MYH4,ACTA1,TNNT3,MYBPC2,TTN",
         "geneSymbolmore2": "MYH11,TAGLN,PECAM1,CD68,COL1A1"},

        {"tissueType": "Muscle", "cellName": "Satellite cells", "shortName": "SatCell",
         "geneSymbolmore1": "PAX7,MYOD1,MYF5,CALCR,VCAM1",
         "geneSymbolmore2": "MYH1,MYH2,PECAM1,CD68"},

        {"tissueType": "Muscle", "cellName": "Smooth muscle cells", "shortName": "SMC",
         "geneSymbolmore1": "ACTA2,MYH11,TAGLN,CNN1,PDGFRB",
         "geneSymbolmore2": "MYH1,MYH2,PECAM1,CD68"},
    ]

    all_cells.extend(muscle_cells)

    #============================================================================#
    # SKIN
    #============================================================================#

    skin_cells = [
        {"tissueType": "Skin", "cellName": "Keratinocytes", "shortName": "Keratino",
         "geneSymbolmore1": "KRT1,KRT10,KRT14,KRT5,DSG3,DSC3",
         "geneSymbolmore2": "PECAM1,CD68,COL1A1,MLANA"},

        {"tissueType": "Skin", "cellName": "Melanocytes", "shortName": "Melano",
         "geneSymbolmore1": "MLANA,PMEL,DCT,TYRP1,TYR,SOX10",
         "geneSymbolmore2": "KRT14,COL1A1,PECAM1,CD68"},

        {"tissueType": "Skin", "cellName": "Fibroblasts", "shortName": "Fibro",
         "geneSymbolmore1": "COL1A1,COL1A2,COL3A1,DCN,LUM,PDGFRA",
         "geneSymbolmore2": "KRT14,PECAM1,CD68,MLANA"},

        {"tissueType": "Skin", "cellName": "Langerhans cells", "shortName": "LC",
         "geneSymbolmore1": "CD207,CD1A,CD1C,EPCAM,ITGAX",
         "geneSymbolmore2": "KRT14,COL1A1,PECAM1,MLANA"},
    ]

    all_cells.extend(skin_cells)

    #============================================================================#
    # ADIPOSE
    #============================================================================#

    adipose_cells = [
        {"tissueType": "Adipose", "cellName": "Adipocytes", "shortName": "Adipo",
         "geneSymbolmore1": "ADIPOQ,LEP,PLIN1,FABP4,LPL,PPARG",
         "geneSymbolmore2": "PECAM1,CD68,COL1A1,ACTA2"},

        {"tissueType": "Adipose", "cellName": "Preadipocytes", "shortName": "PreAdipo",
         "geneSymbolmore1": "PDGFRA,CD34,DLK1,PPARG,C3",
         "geneSymbolmore2": "ADIPOQ,LEP,PLIN1,PECAM1,CD68"},
    ]

    all_cells.extend(adipose_cells)

    # Create DataFrame
    df = pd.DataFrame(all_cells)

    return df

def main():
    print("Creating enhanced ScType database...")

    # Create the database
    enhanced_db = create_enhanced_database()

    # Print summary statistics
    print(f"\nDatabase created with {len(enhanced_db)} entries across {enhanced_db['tissueType'].nunique()} tissue types")

    # Save to Excel file
    output_file = Path("/home/user/sc-type/ScTypeDB_enhanced.xlsx")
    enhanced_db.to_excel(output_file, index=False, engine='openpyxl')

    print(f"Database saved to: {output_file}")

    # Print summary by tissue
    print("\n=== Database Summary ===")
    tissue_summary = enhanced_db.groupby('tissueType').size().sort_index()
    for tissue, count in tissue_summary.items():
        print(f"{tissue:20s}: {count:3d} cell types")

    print("\nDone!")

if __name__ == "__main__":
    main()
