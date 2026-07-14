"""
crispr_designer/sacas9_scores.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The published SaCas9 specificity model coefficients from:

    Tycko, J., Barrera, L.A., Huston, N.C. et al. "Pairwise library screen
    systematically interrogates Staphylococcus aureus Cas9 specificity in human
    cells." Nat Commun 9, 2962 (2018). https://doi.org/10.1038/s41467-018-05391-2

These exact values were retrieved from the paper's official code repository
(github.com/editasmedicine/pairwise-library-screen, models/rstan_mult_21_coeffs.txt) -
fetched directly, not reconstructed from memory, so they can be trusted to match the
original table. Only the 21nt-guide table is included here because
crispr_designer.systems.SaCas9System always uses a fixed 21nt spacer; the source repo
also has separate tables for 19/20/22/23/24nt guides (rstan_mult_<len>_coeffs.txt) since
the model was fit independently per guide length.

Model background: this is a nonlinear multiplicative-penalty model (Bayesian-fit via
Hamiltonian Monte Carlo / Rstan) rather than a simple frequency table like SpCas9's CFD
score, but it plays the same role - a position- and mismatch-identity-specific score in
[0, 1] that multiplies together across every mismatch position to estimate relative
off-target activity (1.0 = same as a perfect match, lower = more disrupted).

SACAS9_MISMATCH_SCORES is keyed '<hit_base><guide_base>,<position>' where position is
1-indexed counting from the spacer's 5' end (PAM-distal, since SaCas9's PAM is 3' like
SpCas9's) - the same direction/convention used by CFD_MISMATCH_SCORES in cfd_scores.py.
This is the reverse of the source file's own "Position" column (which counts from the
PAM-proximal end); the mapping used here is position = 22 - file_position, replicating
predict_activity_single.py's own `par_dict[guide_len][guide_len - Position + 1]` indexing.

IMPORTANT: per the original reference implementation (predict_activity_single.py), position
1 (the very 5'-most spacer base, i.e. the source file's Position=21 rows) is deliberately
never scored - "mismatches at the 5' position are not truly determined by the model."
crispr_designer.offtarget.calc_sacas9_score replicates that by always skipping position 1,
so position-1 entries are omitted below entirely (they would never be looked up).

This table is only valid for a 21nt SaCas9 spacer with an NNGRRT PAM, and only covers
single mismatches (the source model has no bulge/gap support) - see
crispr_designer/offtarget.py's offtarget_risk_score() for how this is enforced.
"""

SACAS9_MISMATCH_SCORES = {
    # position 2 (source file Position 20)
    'AC,2': 0.546734106002513, 'AG,2': 0.992970435672823, 'AT,2': 0.657829492961583,
    'CA,2': 0.522967414281177, 'CG,2': 0.643567704594203, 'CT,2': 0.618178833577651,
    'GA,2': 0.46259280684568, 'GC,2': 0.15914318618166, 'GT,2': 0.490930802795924,
    'TA,2': 0.624360922159452, 'TC,2': 0.640457481900362, 'TG,2': 0.703936911913447,

    # position 3 (source file Position 19)
    'AC,3': 0.547691877828188, 'AG,3': 0.732601715905163, 'AT,3': 0.636416873002814,
    'CA,3': 0.64726672914593, 'CG,3': 0.681666229047031, 'CT,3': 0.581015956193127,
    'GA,3': 0.643038591122457, 'GC,3': 0.274409757488698, 'GT,3': 0.366150894868735,
    'TA,3': 0.698428615250848, 'TC,3': 0.60640601604755, 'TG,3': 0.475908917669472,

    # position 4 (source file Position 18)
    'AC,4': 0.552964194291369, 'AG,4': 0.980427129883399, 'AT,4': 0.363267257971301,
    'CA,4': 0.178476371392183, 'CG,4': 0.416784632985585, 'CT,4': 0.48725485300852,
    'GA,4': 0.372212865767711, 'GC,4': 0.289442475421799, 'GT,4': 0.518288493112618,
    'TA,4': 0.277678942354815, 'TC,4': 0.629270731339457, 'TG,4': 0.168587534931865,

    # position 5 (source file Position 17)
    'AC,5': 0.538528380008684, 'AG,5': 0.685646931401168, 'AT,5': 0.576799741012485,
    'CA,5': 0.589826846023653, 'CG,5': 0.495035619290045, 'CT,5': 0.942657171400989,
    'GA,5': 0.680737138305581, 'GC,5': 0.300672177947638, 'GT,5': 0.796643503923613,
    'TA,5': 0.38993398748401, 'TC,5': 0.390709495232679, 'TG,5': 0.515932076149009,

    # position 6 (source file Position 16)
    'AC,6': 0.325610215099261, 'AG,6': 0.844922889454443, 'AT,6': 0.642009817365257,
    'CA,6': 0.445331765504143, 'CG,6': 0.531756470458885, 'CT,6': 0.540642115587474,
    'GA,6': 0.552943486205603, 'GC,6': 0.368351594652677, 'GT,6': 0.530248054175443,
    'TA,6': 0.422521983677337, 'TC,6': 0.280235604556182, 'TG,6': 0.503183576570521,

    # position 7 (source file Position 15)
    'AC,7': 0.576468070360727, 'AG,7': 0.667118034664953, 'AT,7': 0.567712358758302,
    'CA,7': 0.685794931847784, 'CG,7': 0.676004709924651, 'CT,7': 0.873378872726365,
    'GA,7': 0.584975615068689, 'GC,7': 0.349878390379607, 'GT,7': 0.645199638303918,
    'TA,7': 0.557841348776102, 'TC,7': 0.717284514586225, 'TG,7': 0.565247609715247,

    # position 8 (source file Position 14)
    'AC,8': 0.127378670053041, 'AG,8': 0.683173864400938, 'AT,8': 0.599538413152022,
    'CA,8': 0.325409595348562, 'CG,8': 0.458277575190559, 'CT,8': 0.678015697554753,
    'GA,8': 0.582110823409603, 'GC,8': 0.103481082335252, 'GT,8': 0.51997749852277,
    'TA,8': 0.104010305290341, 'TC,8': 0.41657664122309, 'TG,8': 0.48471892642459,

    # position 9 (source file Position 13)
    'AC,9': 0.633439692869058, 'AG,9': 0.555513095453494, 'AT,9': 0.506919342187349,
    'CA,9': 0.64243781038661, 'CG,9': 0.575785506970191, 'CT,9': 0.58682811500121,
    'GA,9': 0.495378197462189, 'GC,9': 0.529851864152853, 'GT,9': 0.326200207604272,
    'TA,9': 0.441491391101617, 'TC,9': 0.576259924027269, 'TG,9': 0.56059950527931,

    # position 10 (source file Position 12)
    'AC,10': 0.615658340132625, 'AG,10': 0.520509986734207, 'AT,10': 0.646750045500538,
    'CA,10': 0.224012343553474, 'CG,10': 0.444440796592595, 'CT,10': 0.686634819793447,
    'GA,10': 0.27782470178935, 'GC,10': 0.43201684115079, 'GT,10': 0.406797059100196,
    'TA,10': 0.388385789116283, 'TC,10': 0.646131400505293, 'TG,10': 0.46065379532187,

    # position 11 (source file Position 11)
    'AC,11': 0.487865019225753, 'AG,11': 0.457991759790709, 'AT,11': 0.455951859040721,
    'CA,11': 0.217821657175021, 'CG,11': 0.387918368049248, 'CT,11': 0.431002446808513,
    'GA,11': 0.307880431933795, 'GC,11': 0.155692052708749, 'GT,11': 0.189198509598849,
    'TA,11': 0.507615988748957, 'TC,11': 0.555381735098051, 'TG,11': 0.374430944727598,

    # position 12 (source file Position 10)
    'AC,12': 0.42689462494337, 'AG,12': 0.612354705406714, 'AT,12': 0.456786357841618,
    'CA,12': 0.0642029910217309, 'CG,12': 0.216680572369354, 'CT,12': 0.261650340816763,
    'GA,12': 0.299394927832334, 'GC,12': 0.189216008110395, 'GT,12': 0.23819640621691,
    'TA,12': 0.0630783302069826, 'TC,12': 0.251895624458276, 'TG,12': 0.23572096660813,

    # position 13 (source file Position 9)
    'AC,13': 0.102397715418986, 'AG,13': 0.37824405562458, 'AT,13': 0.267574282087781,
    'CA,13': 0.0615709293577342, 'CG,13': 0.228139311830686, 'CT,13': 0.317651503071752,
    'GA,13': 0.141962480894423, 'GC,13': 0.0529315895573772, 'GT,13': 0.0970060141425968,
    'TA,13': 0.0388099660522349, 'TC,13': 0.101076231082529, 'TG,13': 0.0794278765493476,

    # position 14 (source file Position 8)
    'AC,14': 0.200257694187704, 'AG,14': 0.395337246155066, 'AT,14': 0.12373852758353,
    'CA,14': 0.0162972303713068, 'CG,14': 0.0287977020564966, 'CT,14': 0.254249722211379,
    'GA,14': 0.0292199017744031, 'GC,14': 0.039536033512892, 'GT,14': 0.0302548103084032,
    'TA,14': 0.0460124390086707, 'TC,14': 0.346568359640956, 'TG,14': 0.0553640525834663,

    # position 15 (source file Position 7)
    'AC,15': 0.136499862475113, 'AG,15': 0.27104610293669, 'AT,15': 0.0789972001603137,
    'CA,15': 0.0226050940568139, 'CG,15': 0.128405660784272, 'CT,15': 0.187839704652641,
    'GA,15': 0.0195124051391844, 'GC,15': 0.0232766388217564, 'GT,15': 0.0307495242144697,
    'TA,15': 0.0961446136331311, 'TC,15': 0.178707425151124, 'TG,15': 0.0409451753782538,

    # position 16 (source file Position 6)
    'AC,16': 0.0493188093803229, 'AG,16': 0.451540654839527, 'AT,16': 0.10911606478632,
    'CA,16': 0.0597797023253567, 'CG,16': 0.191384782904834, 'CT,16': 0.0979027817488908,
    'GA,16': 0.140212678240853, 'GC,16': 0.0191474738800195, 'GT,16': 0.0522509311690574,
    'TA,16': 0.111029408600699, 'TC,16': 0.117130325394167, 'TG,16': 0.222056916803022,

    # position 17 (source file Position 5)
    'AC,17': 0.136539686900136, 'AG,17': 0.443502321156338, 'AT,17': 0.199571974143948,
    'CA,17': 0.050358810530371, 'CG,17': 0.0357634305414073, 'CT,17': 0.236699959495542,
    'GA,17': 0.107550670622112, 'GC,17': 0.025852157724115, 'GT,17': 0.0575108872220017,
    'TA,17': 0.056782864029314, 'TC,17': 0.117154009003686, 'TG,17': 0.0586170001847033,

    # position 18 (source file Position 4)
    'AC,18': 0.0927354706350068, 'AG,18': 0.584714691584839, 'AT,18': 0.16959064954525,
    'CA,18': 0.0775368431850832, 'CG,18': 0.0367789926238243, 'CT,18': 0.113816114093481,
    'GA,18': 0.0533316085854886, 'GC,18': 0.0148089478006535, 'GT,18': 0.0400655724089505,
    'TA,18': 0.093036719073791, 'TC,18': 0.0979767483566834, 'TG,18': 0.0832194891199737,

    # position 19 (source file Position 3)
    'AC,19': 0.130261909041357, 'AG,19': 0.436629040331668, 'AT,19': 0.296554470273688,
    'CA,19': 0.0152863394236489, 'CG,19': 0.070695309584881, 'CT,19': 0.0573652386585105,
    'GA,19': 0.075377597399379, 'GC,19': 0.0356756456096336, 'GT,19': 0.0554196638437215,
    'TA,19': 0.059627368313025, 'TC,19': 0.136926154119441, 'TG,19': 0.0372902912891327,

    # position 20 (source file Position 2)
    'AC,20': 0.0773508496492739, 'AG,20': 0.22354371179448, 'AT,20': 0.0601305506668423,
    'CA,20': 0.1432549909342, 'CG,20': 0.158363763886502, 'CT,20': 0.157638078782083,
    'GA,20': 0.0460670515595503, 'GC,20': 0.0154733404301811, 'GT,20': 0.0535733539298569,
    'TA,20': 0.218904216245295, 'TC,20': 0.191997864741009, 'TG,20': 0.122380573189007,

    # position 21 (source file Position 1)
    'AC,21': 0.105693733152788, 'AG,21': 0.445155052527266, 'AT,21': 0.0862142526234551,
    'CA,21': 0.0426130973647857, 'CG,21': 0.104525045604611, 'CT,21': 0.179929057922201,
    'GA,21': 0.109313834975141, 'GC,21': 0.0449961433951072, 'GT,21': 0.0805457311950328,
    'TA,21': 0.0638608590837734, 'TC,21': 0.043097315653776, 'TG,21': 0.177723354855625,
}
