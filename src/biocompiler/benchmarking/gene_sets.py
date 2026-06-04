"""
BioCompiler Benchmarking Gene Sets
====================================
Curated gene sets for benchmarking protein optimization pipelines.

Each gene set maps a gene symbol to a dict with:
  - protein_sequence : str   – amino acid sequence (UniProt canonical, with signal peptide)
  - organism         : str   – source organism
  - description      : str   – human-readable protein name
  - uniprot_id       : str   – UniProtKB accession (or empty string if unknown)

Sequences sourced from UniProtKB/Swiss-Prot (reviewed, canonical isoforms).
"""

from __future__ import annotations

__all__ = [
    "HUMAN_THERAPEUTIC_GENES",
    "VACCINE_ANTIGEN_GENES",
    "STRESS_TEST_GENES",
    "get_all_gene_sets",
]


# except where noted (e.g. F8 truncated to first 300 aa).

HUMAN_THERAPEUTIC_GENES: dict[str, dict] = {
    "INS": {
        "protein_sequence": (
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
            "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLEN"
            "YCN"
        ),
        "organism": "Homo sapiens",
        "description": "Insulin precursor (preproinsulin); therapeutic: recombinant human insulin",
        "uniprot_id": "P01308",
    },
    "GH1": {
        "protein_sequence": (
            "MATGSRTSLLLAFGLLCLPWLQEGSAFPTIPLSRLFDNAMLRAHRLHQLAFDT"
            "YQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISL"
            "LLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSP"
            "RTGQIFKQTYSKFDTNSHNDDALLKNYGLLYCFRKDMDKVETFLRIVQCRSVE"
            "GSCGF"
        ),
        "organism": "Homo sapiens",
        "description": "Somatotropin (growth hormone); therapeutic: recombinant human GH (somatropin)",
        "uniprot_id": "P01241",
    },
    "EPO": {
        "protein_sequence": (
            "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITT"
            "GCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQAL"
            "LVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTI"
            "TADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR"
        ),
        "organism": "Homo sapiens",
        "description": "Erythropoietin precursor; therapeutic: recombinant epoetin alfa",
        "uniprot_id": "P01588",
    },
    "IFNA2": {
        "protein_sequence": (
            "MALTFALLVALLVLSCKSSCSVGCDLPQTHSLGSRRTLMLLAQMRRISLFSCLK"
            "DRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTKDSSAAWDETLLDKF"
            "YTELYQQLNDLEACVIQGVGVTETPLMKEDSILAVRKYFQRITLYLKEKKYSP"
            "CAWEVVRAEIMRSFSLSTNLQESLRSKE"
        ),
        "organism": "Homo sapiens",
        "description": "Interferon alpha-2 precursor; therapeutic: recombinant IFN-alpha-2a/2b",
        "uniprot_id": "P01563",
    },
    "IL2": {
        "protein_sequence": (
            "MYRMQLLSCIALSLALVTNSAPTSSSTKKTQLQLEHLLLDLQMILNGINNYKNP"
            "KLTRMLTFKFYMPKKATELKHLQCLEEELKPLEEVLNLAQSKNFHLRPRDLISN"
            "INVIVLELKGSETTFMCEYADETATIVEFLNRWITFCQSIISTLT"
        ),
        "organism": "Homo sapiens",
        "description": "Interleukin-2 precursor; therapeutic: recombinant aldesleukin",
        "uniprot_id": "P60568",
    },
    "CSF3": {
        "protein_sequence": (
            "MAGPATQSPMKLMALQLLLWHSALWTVQEATPLGPASSLPQSFLLKCLEQVRKI"
            "QGDGAALQEKLVSECATYKLCHPEELVLLGHSLGIPWAPLSSCPSQALQLAGCL"
            "SQLHSGLFLYQGLLQALEGISPELGPTLDTLQLDVADFATTIWQQMEELGMAPA"
            "LQPTQGAMPAFASAFQRRAGGVLVASHLQSFLEVSYRVLRHLAQP"
        ),
        "organism": "Homo sapiens",
        "description": "Granulocyte colony-stimulating factor precursor; therapeutic: filgrastim",
        "uniprot_id": "P09919",
    },
    "PLAT": {
        "protein_sequence": (
            "MDAMKRGLCCVLLLCGAVFVSPSQEIHARFRRGARSYQVICRDEKTQMIYQQHQ"
            "SWLRPVLRSNRVEYCWCNSGRAQCHSVPVKSCSEPRCFNGGTCQQALYFSDFVC"
            "QCPEGFAGKCCEIDTRATCYEDQGISYRGTWSTAESGAECTNWNSSALAQKPYS"
            "GRRPDAIRLGLGNHNYCRNPDRDSKPWCYVFKAGKYSSEFCSTPACSEGNSDCY"
            "FGNGSAYRGTHSLTESGASCLPWNSMILIGKVYTAQNPSAQALGLGKHNYCRNP"
            "DGDAKPWCHVLKNRRLTWEYCDVPSCSTCGLRQYSQPQFRIKGGLFADIASHPW"
            "QAAIFAKHRRSPGERFLCGGILISSCWILSAAHCFQERFPPHHLTVILGRTYRV"
            "VPGEEEQKFEVEKYIVHKEFDDDTYDNDIALLQLKSDSSRCAQESSVVRTVCLP"
            "PADLQLPDWTECELSGYGKHEALSPFYSERLKEAHVRLYPSSRCTSQHLLNRTV"
            "TDNMLCAGDTRSGGPQANLHDACQGDSGGPLVCLNDGRMTLVGIISWGLGCGQK"
            "DVPGVYTKVTNYLDWIRDNMRP"
        ),
        "organism": "Homo sapiens",
        "description": "Tissue-type plasminogen activator precursor; therapeutic: alteplase",
        "uniprot_id": "P00750",
    },
    "F8": {
        # First 300 aa of the 2351-aa precursor (A1 domain region)
        "protein_sequence": (
            "MQIELSTCFFLCLLRFCFSATRRYYLGAVELSWDYMQSDLGELPVDARFPPRVPK"
            "SFPFNTSVVYKKTLFVEFTDHLFNIAKPRPPWMGLLGPTIQAEVYDTVVITLKN"
            "MASHPVSLHAVGVSYWKASEGAEYDDQTSQREKEDDKVFPGGSHTYVWQVLKENG"
            "PMASDPLCLTYSYLSHVDLVKDLNSGLIGALLVCREGSLAKEKTQTLHKFILLFA"
            "VFDEGKSWHSETKNSLMQDRDAASARAWPKMHTVNGYVNRSLPGLIGCHRKSVYW"
            "HVIGMGTTPEVHSIFLEGHTFLVRNH"
        ),
        "organism": "Homo sapiens",
        "description": "Coagulation factor VIII precursor (first 300 aa); therapeutic: recombinant FVIII",
        "uniprot_id": "P00451",
    },
    "SERPINA1": {
        "protein_sequence": (
            "MPSSVSWGILLLAGLCCLVPVSLAEDPQGDAAQKTDTSHHDQDHPTFNKITPNLA"
            "EFAFSLYRQLAHQSNSTNIFFSPVSIATAFAMLSLGTKADTHDEILEGLNFNLTE"
            "IPEAQIHEGFQELLRTLNQPDSQLQLTTGNGLFLSEGLKLVDKFLEDVKKLYHSE"
            "AFTVNFGDTEEAKKQINDYVEKGTQGKIVDLVKELDRDTVFALVNYIFFKGKWER"
            "PFEVKDTEEEDFHVDQVTTVKVPMMKRLGMFNIQHCKKLSSWVLLMKYLGNATAI"
            "FFLPDEGKLQHLENELTHDIITKFLENEDRRSASLHLPKLSITGTYDLKSVLGQL"
            "GITKVFSNGADLSGVTEEAPLKLSKAVHKAVLTIDEKGTEAAGAMFLEAIPMSIPP"
            "EVKFNKPFVFLMIEQNTKSPLFMGKVVNPTQK"
        ),
        "organism": "Homo sapiens",
        "description": "Alpha-1-antitrypsin precursor; therapeutic: recombinant A1AT augmentation therapy",
        "uniprot_id": "P01009",
    },
    "GBA": {
        "protein_sequence": (
            "MEFSSPSREECPKPLSRVSIMAGSLTGLLLLQAVSWASGARPCIPKSFGYSSVVC"
            "VCNATYCDSFDPPTFPALGTFSRYESTRSGRRMELSMGPIQANHTGTGLLLTLQP"
            "EQKFQKVKGFGGAMTDAAALNILALSPPAQNLLLKSYFSEEGIGYNIIRVPMASC"
            "DFSIRTYTYADTPDDFQLHNFSLPEEDTKLKIPLIHRALQLAQRPVSLLASPWTSP"
            "TWLKTNGAVNGKGSLKGQPGDIYHQTWARYFVKFLDAYAEHKLQFWAVTAENEPS"
            "AGLLSGYPFQCLGFTPEHQRDFIARDLGPTLANSTHHNVRLLMLDDQRLLLPHWA"
            "KVVLTDPEAAKYVHGIAVHWYLDFLAPAKATLGETHRLFPNTMLFASEACVGSKF"
            "WEQSVRLGSWDRGMQYSHSIITNLLYHVVGWTDWNLALNPEGGPNWVRNFVDSP"
            "IIVDITKDTFYKQPMFYHLGHFSKFIPEGSQRVGLVASQKNDLDAVALMHPDGSA"
            "VVVVLNRSSKDVPLTIKDPAVGFLETISPGYSIHTYLWRRQ"
        ),
        "organism": "Homo sapiens",
        "description": "Lysosomal acid glucosylceramidase precursor; therapeutic: imiglucerase/velaglucerase",
        "uniprot_id": "P04062",
    },
}


# ---------------------------------------------------------------------------
# Vaccine Antigen Genes
# ---------------------------------------------------------------------------
# Pathogen-derived antigens used in licensed or experimental vaccines.
# Full UniProt canonical sequences where available.
# Includes all 10 task-specified vaccine antigens plus additional entries.

VACCINE_ANTIGEN_GENES: dict[str, dict] = {
    "SARS2_Spike": {
        "protein_sequence": (
            "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS"
            "NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV"
            "NNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLE"
            "GKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQT"
            "LLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETK"
            "CTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISN"
            "CVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIAD"
            "YNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPC"
            "NGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN"
            "FNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITP"
            "GTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSY"
            "ECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTI"
            "SVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQE"
            "VFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDC"
            "LGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAM"
            "QMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALN"
            "TLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRA"
            "SANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPA"
            "ICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDP"
            "LQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDL"
            "QELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDD"
            "SEPVLKGVKLHYT"
        ),
        "organism": "Severe acute respiratory syndrome coronavirus 2",
        "description": "SARS-CoV-2 Spike glycoprotein (Wuhan-Hu-1); vaccine: BNT162b2/mRNA-1273",
        "uniprot_id": "P0DTC2",
    },

    "SARS2_RBD": {
        "protein_sequence": (
            "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFK"
            "CYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNS"
            "NNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQ"
            "PTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
        ),
        "organism": "Severe acute respiratory syndrome coronavirus 2",
        "description": "SARS-CoV-2 Spike RBD domain (residues 319-541); vaccine: subunit vaccines",
        "uniprot_id": "P0DTC2",
    },

    "H1N1_HA": {
        "protein_sequence": (
            "MKANLLVLLCALAAADADTICIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDSHNGKLCR"
            "LKGIAPLQLGKCNIAGWLLGNPECDPLLPVRSWSYIVETPNSENGICYPGDFIDYEELRE"
            "QLSSVSSFERFEIFPKESSWPNHNTNGVTAACSHEGKSSFYRNLLWLTEKEGSYPKLKNS"
            "YVNKKGKEVLVLWGIHHPPNSKEQQNLYQNENAYVSVVTSNYNRRFTPEIAERPKVRDQA"
            "GRMNYYWTLLKPGDTIIFEANGNLIAPMYAFALSRGFGSGIITSNASMHECNTKCQTPLG"
            "AINSSLPYQNIHPVTIGECPKYVRSAKLRMVTGLRNIPSIQSRGLFGAIAGFIEGGWTGM"
            "IDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNTVIEKMNIQFTAVGKEFNKLEKRM"
            "ENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDSNVKNLYEKVKSQLKNNAKEIGNGC"
            "FEFYHKCDNECMESVRNGTYDYPKYSEESKLNREKVDGVKLESMGIYQILAIYSTVASSL"
            "VLLVSLGAISFWMCSNGSLQCRICI"
        ),
        "organism": "Influenza A virus (A/Puerto Rico/8/1934 H1N1)",
        "description": "Influenza A H1N1 Hemagglutinin; vaccine: seasonal flu vaccines",
        "uniprot_id": "P03452",
    },

    "H3N2_HA": {
        "protein_sequence": (
            "MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQ"
            "SSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPD"
            "YASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVL"
            "NVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVR"
            "GLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITP"
            "NGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWE"
            "GMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEG"
            "RIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGN"
            "GCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISC"
            "FLLCVVLLGFIMWACQRGNIRCNICI"
        ),
        "organism": "Influenza A virus (A/Aichi/2/1968 H3N2)",
        "description": "Influenza A H3N2 Hemagglutinin; vaccine: seasonal flu vaccines",
        "uniprot_id": "P03437",
    },

    "HIV1_gp120": {
        "protein_sequence": (
            "TEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVT"
            "ENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIME"
            "KGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSF"
            "EPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEV"
            "VIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAH"
            "CNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQL"
            "FNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNIT"
            "GLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREK"
            "RAVGIGALFLGFLGAAGSTM"
        ),
        "organism": "Human immunodeficiency virus type 1",
        "description": "HIV-1 Envelope gp120 (residues 31-530 of gp160); vaccine: experimental HIV vaccines",
        "uniprot_id": "P04578",
    },

    "RSV_F": {
        "protein_sequence": (
            "MELLILKANAITTILTAVTFCFASGQNITEEFYQSTCSAVSKGYLSALRTGWYTSVITIE"
            "LSNIKENKCNGTDAKVKLIKQELDKYKNAVTELQLLMQSTPPTNNRARRELPRFMNYTLN"
            "NAKKTNVTLSKKRKRRFLGFLLGVGSAIASGVAVSKVLHLEGEVNKIKSALLSTNKAVVS"
            "LSNGVSVLTSKVLDLKNYIDKQLLPIVNKQSCSISNIETVIEFQQKNNRLLEITREFSVN"
            "AGVTTPVSTYMLTNSELLSLINDMPITNDQKKLMSNNVQIVRQQSYSIMSIIKEEVLAYV"
            "VQLPLYGVIDTPCWKLHTSPLCTTNTKEGSNICLTRTDRGWYCDNAGSVSFFPQAETCKV"
            "QSNRVFCDTMNSLTLPSEINLCNVDIFNPKYDCKIMTSKTDVSSSVITSLGAIVSCYGKT"
            "KCTASNKNRGIIKTFSNGCDYVSNKGMDTVSVGNTLYYVNKQEGKSLYVKGEPIINFYDP"
            "LVFPSDEFDASISQVNEKINQSLAFIRKSDELLHNVNAGKSTTNIMITTIIIVIIVILLS"
            "LIAVGLLLYCKARSTPVTLSKDQLSGINNIAFSN"
        ),
        "organism": "Human respiratory syncytial virus",
        "description": "RSV Fusion glycoprotein F0; vaccine: Arexvy/Abrysvo",
        "uniprot_id": "P03420",
    },

    "RABV_G": {
        "protein_sequence": (
            "MSKIFVNPSAIRAGLADLEMAEETVDLINRNIEDNDAHLQGEPIEVDNLPEDMKRLHLDD"
            "EKSSNLGEMVRVGEGKYREDFQMDEGEDPNLLFQSYLDNVGVQIVRQMRSGERFLKIWSQ"
            "TVEEIVSYVTVNFPNPPRRSSEDKSTQTTGRELKKETTSAFSQRESQPSKARMVAQVAPG"
            "PPALEWSATNEEDDLSVEAEIAHQIAESFSKKYKFPSRSSGIFLYNFEQLKMNLDDIVKE"
            "AKNVPGVTRLAHDGSKIPLRCVLGWVALANSKKFQLLVEADKLSKIMQDDLNRYTSC"
        ),
        "organism": "Rabies virus (strain PV)",
        "description": "Glycoprotein G; vaccine: HDCV/PCECV rabies vaccines",
        "uniprot_id": "P15198",
    },

    "Zika_E": {
        "protein_sequence": (
            "IRCIGVSNRDFVEGMSGGTWVDVVLEHGGCVTVMAQDKPTVDIELVTTTVSNMAEVRSYC"
            "YEASISDMASDSRCPTQGEAYLDKQSDTQYVCKRTLVDRGWGNGCGLFGKGSLVTCAKFT"
            "CSKKMTGKSIQPENLEYRIMLSVHGSQHSGMIGYETDEDRAKVEVTPNSPRAEATLGGFG"
            "SLGLDCEPRTGLDFSDLYYLTMNNKHWLVHKEWFHDIPLPWHAGADTGTPHWNNKEALVE"
            "FKDAHAKRQTVVVLGSQEGAVHTALAGALEAEMDGAKGRLFSGHLKCRLKMDKLRLKGVS"
            "YSLCTAAFTFTKVPAETLHGTVTVEVQYAGTDGPCKIPVQMAVDMQTLTPVGRLITANPV"
            "ITESTENSKMMLELDPPFGDSYIVIGVGDKKITHHWHRSGSTIGKAFEATVRGAKRMAVL"
            "GDTAWDFGSVGGVFNSLGKGIHQIFGAAFKSLFGGMSWFSQILIGTLLVWLGLNTKNGSI"
            "SLTCLALGGVMIFLSTAVSA"
        ),
        "organism": "Zika virus",
        "description": "Zika virus Envelope protein E (residues 291-790 of polyprotein); vaccine: experimental",
        "uniprot_id": "Q32ZE1",
    },

    "DENV2_E": {
        "protein_sequence": (
            "MRCIGISNRDFVEGVSGGSWVDIVLEHGSCVTTMAKNKPTLDFELIETEAKQPATLRKYC"
            "IEAKLTNTTTDSRCPTQGEPSLNEEQDKRFVCKHSMVDRGWGNGCGLFGKGGIVTCAMFT"
            "CKKNMKGKVVQPENLEYTIVITPHSGEEHAVGNDTGKHGKEIKITPQSSITEAELTGYGT"
            "VTMECSPRTGLDFNEMVLLQMENKAWLVHRQWFLDLPLPWLPGADTQGSNWIQKETLVTF"
            "KNPHAKKQDVVVLGSQEGAMHTALTGATEIQMSSGNLLFTGHLKCRLRMDKLQLKGMSYS"
            "MCTGKFKVVKEIAETQHGTIVIRVQYEGDGSPCKIPFEIMDLEKRHVLGRLITVNPIVTE"
            "KDSPVNIEAEPPFGDSYIIIGVEPGQLKLNWFKKGSSIGQMIETTMRGAKRMAILGDTAW"
            "DFGSLGGVFTSIGKALHQVFGAIYGAAFSGVSWIMKILIGVIITWIGMNSRSTSLSVSLV"
            "LVGVVTLYLGVMVQA"
        ),
        "organism": "Dengue virus type 2",
        "description": "Dengue virus type 2 Envelope protein E (residues 281-775 of polyprotein); vaccine: Dengvaxia",
        "uniprot_id": "P14340",
    },

    "MTB_Ag85B": {
        "protein_sequence": (
            "MPRVEVGLVIHSRMHARAPVDVWRSVRSLPDFWRLLQVRVASQFGDGLFQAGLAGALLFN"
            "PDRAADPMAIAGAFAVLFLPYSLLGPFAGALMDRWDRRWVLVGANTGRLALIAGVGTILA"
            "VGAGDVPLLVGALVANGLARFVASGLSAALPHVVPREQVVTMNSVAIASGAVSAFLGANF"
            "MLLPRWLLGSGDEGASAIVFLVAIPVSIALLWSLRFGPRVLGPDDTERAIHGSAVYAVVT"
            "GWLHGARTVVQLPTVAAGLSGLAAHRMVVGINSLLILLLVRHVTARAVGGLGTALLFFAA"
            "TGLGAFLANVLTPTAIRRWGRYATANGALAAAATIQVAAAGLLVPVMVVCGFLLGVAGQV"
            "VKLCADSAMQMDVDDALRGHVFAVQDALFWVSYILSITVAAALIPEHGHAPVFVLFGSAI"
            "YLAGLVVHTIVGRRGQPVIGR"
        ),
        "organism": "Mycobacterium tuberculosis",
        "description": "Antigen 85B; vaccine: M72/AS01E candidate",
        "uniprot_id": "P0A5C2",
    },

    "PF_CSP": {
        "protein_sequence": (
            "MMRKLAILSVSSFLFVEALFQEYQCYGSSSNTRVLNELNYDNAGTNLYNELEMNYYGKQE"
            "NWYSLKKNSRSLGENDDGNNEDNEKLRKPKHKKLKQPADGNPDPNANPNVDPNANPNVDP"
            "NANPNVDPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANP"
            "NANPNANPNANPNANPNVDPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANP"
            "NANPNANPNANPNANPNANPNANPNANPNANPNKNNQGNGQGHNMPNDPNRNVDENANAN"
            "SAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQVRIKPGSANKPKDELDY"
            "ANDIEKKICKMEKCSSVFNVVNSSIGLIMVLSFLFLN"
        ),
        "organism": "Plasmodium falciparum",
        "description": "Circumsporozoite protein; vaccine: RTS,S/AS01 (Mosquirix)",
        "uniprot_id": "P19597",
    },
}


# ---------------------------------------------------------------------------
# Stress-Test Genes
# ---------------------------------------------------------------------------
# Synthetic sequences designed to stress-test edge cases in optimization:
#   extreme GC/AT composition, homopolymer runs, restriction site density, etc.

STRESS_TEST_GENES: dict[str, dict] = {
    "PolyA_100": {
        "protein_sequence": (
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        ),
        "organism": "synthetic",
        "description": "Poly-alanine repeat protein (100xA); tests homopolymer codon handling",
        "uniprot_id": "",
        "stress_category": "homopolymer",
    },

    "ExtremeGC_200": {
        "protein_sequence": (
            "RGAARRVGPGAPGGRPAPVGVRAGVAGGVRVRPVAPVRVPRGAAGAGARAAAAVRRGRGA"
            "VRPRVRAGAAAVVARAVPAAPAPVAAVPGGGRRPGAVPVVGRRPARGPPVVAPGVVARRG"
            "RPRPGAGVVVAGVVGPGRRGPPAVPAPRAAVRPPGAAPAAGRAVVGARAGVPPRVGGPPP"
            "RRVGPAVAGPPAAVPVPGVV"
        ),
        "organism": "synthetic",
        "description": "Extreme GC-biased protein (200aa, ~80% GC codons); tests high-GC optimization",
        "uniprot_id": "",
        "stress_category": "high_gc",
    },

    "ExtremeAT_200": {
        "protein_sequence": (
            "LLMKFKNINLKMFNLLFMKNLYFMIYFMKLIKYYKILKYYNYFLKMKNYNNMIFLLIKKL"
            "LLNKMYYLYIMNYLIIIFKYMKYKFLMINKIMYKMFNIILLFMYKINLMLLKKIMLNYKL"
            "LLMLIMNLMLMKMILLKMFNMKIKNLLKYNLFLILKLKIFMNYFNYKYILFMMYNIYNNM"
            "MNNNNNINLYLYNFINFYYF"
        ),
        "organism": "synthetic",
        "description": "Extreme AT-biased protein (200aa, ~20% GC codons); tests low-GC optimization",
        "uniprot_id": "",
        "stress_category": "high_at",
    },

    "Long_1000": {
        "protein_sequence": (
            "IDIYSNRLKESPMDDAQITAGKDFICMAPSRHNVLAQTISMVKVWKPLSLWPIQQKPLKQ"
            "AEQELNCFFAGLHDHYLFTPGAGPTVHPLKDDHAEFGNIRKAGLAKRIAGDEPGPSQAAL"
            "VGKSKPVDANRMLDFPNKPRSVRAEFVGFAFKRSTKRVSGHGLETSRVKSFHPVQSLIHP"
            "QTSRTCIKGATDDSETFQAVRPIPFKISRDFELILHPEFQQRDRLPDDIPVLERDDDALQ"
            "VLFRTITLIIDAVISGAKANLEVRPKSDLIITVEKKPVDDMNNAAILECFFIFYSKMVHL"
            "KYSANEPLERFFVLIAAHEQRVTLVVMQAPIKYDLKIRGEMTEFEVSLIEEVGLVMLGDG"
            "PKPDLYIFRVKCNLLRVYIRSLFLQSVCLGYLTTQSGVLSEELVINLGLLIEDFLCSAPF"
            "SGMSVIALLTTILIPGMNIVEFKTTTGEPPTAAKYVPHRKTIFYASEGNLLSNATKCYNI"
            "QQASYFNIGLKHQTDGICLLVEVMQLFPAKLYESKGKKFPMTKYKEILPREPNRLVLYSQ"
            "SISPERVNVTPVVSGSALIMMLSVRCATLYDASYEGNVEIRTRIMPVHPNGRIGVAGIVT"
            "RFLVGGDRTIMKKASCRLAPDDTENYATNRQNMAASCENAPLPISLRLVASTLIVAIGFI"
            "YAAKAQKLEVLIKRFTKSPDIDEIHRHMGNYASRRPSLLLERGMPTPVKGLDMDRTGPLS"
            "TVSKKAYSQDFEEATLGNKAPDHEGPSLRANKMQMIIQPRHEFIMAGKAEDNQEAPNGFL"
            "SSATAAYTLLFHRAESKSYRTCLLGEEENLRLYGPSSCNDKPMDSIPPIYLKKTKNAIEQ"
            "ALEIASATTKLILSWISKEINARYEFLNQIILNGEGEKSKARFEPMYTELNGVFGKVEGS"
            "SMSPFLKHGGDDVIDNATRSTTAEPRGNSSSNIGMCIEPILESITGVVIEGLVSSNQKTL"
            "RTTEPQYNIVDIGAALCVGAYTKSNESFNFIAAQSLPLRE"
        ),
        "organism": "synthetic",
        "description": "Long synthetic protein (1000aa); tests scalability of optimization pipeline",
        "uniprot_id": "",
        "stress_category": "long_sequence",
    },

    "Repetitive_10x20": {
        "protein_sequence": (
            "LKQFSLEQIAKRFNKMASTALKQFSLEQIAKRFNKMASTALKQFSLEQIAKRFNKMASTA"
            "LKQFSLEQIAKRFNKMASTALKQFSLEQIAKRFNKMASTALKQFSLEQIAKRFNKMASTA"
            "LKQFSLEQIAKRFNKMASTALKQFSLEQIAKRFNKMASTALKQFSLEQIAKRFNKMASTA"
            "LKQFSLEQIAKRFNKMASTA"
        ),
        "organism": "synthetic",
        "description": "Repetitive domain protein (10x 20aa repeat); tests repeat/recombination handling",
        "uniprot_id": "",
        "stress_category": "repetitive",
    },
}


# ---------------------------------------------------------------------------
# Aggregation helper
# ---------------------------------------------------------------------------

def get_all_gene_sets() -> dict[str, dict]:
    """Return all gene sets merged into a single dict.

    In case of key collisions, later sets overwrite earlier ones.
    The merge order is:
      HUMAN_THERAPEUTIC_GENES -> VACCINE_ANTIGEN_GENES -> STRESS_TEST_GENES.
    """
    merged: dict[str, dict] = {}
    merged.update(HUMAN_THERAPEUTIC_GENES)
    merged.update(VACCINE_ANTIGEN_GENES)
    merged.update(STRESS_TEST_GENES)
    return merged


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _validate_sequences() -> None:
    """Validate all gene set sequences at module load time.

    Checks that every sequence:
      1. Contains only valid single-letter amino acid codes
    """
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    errors: list[str] = []

    for gene_set_name, gene_set in [
        ("HUMAN_THERAPEUTIC_GENES", HUMAN_THERAPEUTIC_GENES),
        ("VACCINE_ANTIGEN_GENES", VACCINE_ANTIGEN_GENES),
        ("STRESS_TEST_GENES", STRESS_TEST_GENES),
    ]:
        for name, entry in gene_set.items():
            seq = entry.get("protein_sequence", "")
            if not seq:
                errors.append(
                    f"{gene_set_name}[{name}]: empty or missing sequence"
                )
                continue

            # Check for invalid characters
            invalid = set(seq) - valid_aa
            if invalid:
                errors.append(
                    f"{gene_set_name}[{name}]: invalid amino acids: {invalid}"
                )

    if errors:
        import logging
        logger = logging.getLogger(__name__)
        for err in errors:
            logger.error("Gene set validation error: %s", err)
        raise ValueError(
            f"Gene set validation failed with {len(errors)} error(s): "
            + "; ".join(errors)
        )


# Run validation on import (fail-fast for incorrect data)
_validate_sequences()
