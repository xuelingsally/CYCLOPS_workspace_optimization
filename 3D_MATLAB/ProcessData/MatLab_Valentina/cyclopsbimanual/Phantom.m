% Plot of the Phantom stylus

clear all
% close all
clc

% Phantom geometry
% Number of degrees-of-freedom of the Phantom
n_dof=6;
% Distance between the base frame and the first joint of the Phantom [mm]
d1=168.35;
% Length of the first link of the Phantom [mm]
L1=133.35;
% Length of the second link of the Phantom [mm]
L2=L1;

% THESE PARAMETERS ARE JUST FOR ME TO DO THE PLOTTING
% Radius of the base [mm]
r_base=5;
% Radius of the stylus [mm]
r_stylus=5;
% Radius of the links [mm]
r_link=5;
% Number of slices along the object length
no_slice=10;
% Number of verteces on the object base
no_arc=50;

% Assign the values of the Phantom joint angles
% q=zeros(1,6);

% q=[0.437516, 1.06821, 0.669501, 3.08471, -4.57401, 2.46214;
% 0.437104, 1.06821, 0.669093, 3.08598, -4.57273, 2.46214;
% 0.435871, 1.06821, 0.669093, 3.08471, -4.57017, 2.46214;
% 0.433815, 1.06821, 0.669093, 3.08471, -4.56889, 2.46214;
% 0.430114, 1.06821, 0.669093, 3.08598, -4.56889, 2.46342;
% 0.422301, 1.06821, 0.669093, 3.08471, -4.57017, 2.46342;
% 0.41161, 1.06821, 0.667871, 3.08343, -4.57145, 2.46214;
% 0.404209, 1.06821, 0.666241, 3.07959, -4.57145, 2.46342;
% 0.386938, 1.06821, 0.660128, 3.06809, -4.57145, 2.46342;
% 0.364733, 1.06985, 0.654831, 3.04508, -4.57145, 2.4647;
% 0.337183, 1.0727, 0.650756, 3.02334, -4.57145, 2.48004;
% 0.304287, 1.07637, 0.647496, 2.9965, -4.57145, 2.5056;
% 0.259055, 1.08331, 0.646681, 2.95431, -4.57145, 2.5414;
% 0.210945, 1.08821, 0.645051, 2.90446, -4.57401, 2.5887;
% 0.154611, 1.0931, 0.644644, 2.86355, -4.5804, 2.63088;
% 0.0892302, 1.09555, 0.643014, 2.80091, -4.58551, 2.6833;
% 0.0213823, 1.09555, 0.640976, 2.71142, -4.58807, 2.7421;
% -0.0481103, 1.09555, 0.640976, 2.61427, -4.58679, 2.79452;
% -0.12706, 1.09474, 0.640976, 2.52094, -4.58679, 2.85076;
% -0.197787, 1.09433, 0.643014, 2.43274, -4.58551, 2.91085;
% -0.271391, 1.08984, 0.643014, 2.35731, -4.5804, 2.96326;
% -0.345818, 1.07882, 0.643421, 2.27805, -4.56761, 3.02079;
% -0.420656, 1.05597, 0.649941, 2.1924, -4.54972, 3.08726;
% -0.484803, 1.0168, 0.650756, 2.11314, -4.52287, 3.15246;
% -0.515232, 0.9862, 0.650756, 2.05306, -4.49603, 3.20232;
% -0.515643, 0.984976, 0.650756, 2.03389, -4.48324, 3.22533;
% -0.515643, 0.984976, 0.651164, 2.03133, -4.48324, 3.22788;
% -0.515232, 0.98416, 0.651164, 2.03133, -4.48452, 3.22788;
% -0.51441, 0.974776, 0.650756, 2.03516, -4.4858, 3.23044;
% -0.513176, 0.939277, 0.630382, 2.03389, -4.48708, 3.23811;
% -0.509064, 0.884602, 0.597375, 2.03261, -4.47685, 3.24067;
% -0.504541, 0.818501, 0.551329, 2.03389, -4.45128, 3.24067;
% -0.502896, 0.748729, 0.504061, 2.04156, -4.41549, 3.24067;
% -0.502896, 0.679772, 0.464127, 2.05817, -4.37202, 3.24067;
% -0.502896, 0.607552, 0.427453, 2.06968, -4.33367, 3.24067;
% -0.502896, 0.52717, 0.396892, 2.07863, -4.31066, 3.23939;
% -0.503307, 0.446381, 0.36959, 2.08374, -4.27998, 3.23939;
% -0.503307, 0.372528, 0.343511, 2.08886, -4.26208, 3.23811;
% -0.503307, 0.306428, 0.322322, 2.09397, -4.25569, 3.23939;
% -0.503307, 0.24808, 0.310097, 2.10036, -4.25313, 3.23939;
% -0.503719, 0.186876, 0.302355, 2.11059, -4.25313, 3.24067;
% -0.503719, 0.146482, 0.300725, 2.11442, -4.25313, 3.24067;
% -0.503719, 0.132609, 0.298688, 2.11442, -4.25825, 3.2445;
% -0.503719, 0.131385, 0.294205, 2.11698, -4.26336, 3.24578;
% -0.501251, 0.131385, 0.29013, 2.11826, -4.26208, 3.24578;
% -0.489738, 0.133017, 0.289316, 2.11954, -4.26208, 3.24578;
% -0.463832, 0.137097, 0.281573, 2.12976, -4.26208, 3.2445;
% -0.422301, 0.149746, 0.265681, 2.15917, -4.26208, 3.2445;
% -0.3635, 0.171371, 0.242862, 2.2103, -4.26081, 3.24067;
% -0.28825, 0.184836, 0.21556, 2.28828, -4.25441, 3.21766;
% -0.203543, 0.192997, 0.195186, 2.37138, -4.23652, 3.19337;
% -0.117603, 0.195037, 0.184184, 2.45319, -4.2199, 3.1678;
% -0.0238495, 0.201973, 0.180517, 2.54523, -4.21095, 3.1384;
% 0.072371, 0.208502, 0.181739, 2.65773, -4.21095, 3.11794;
% 0.172292, 0.21095, 0.186221, 2.77918, -4.21095, 3.09493;
% 0.26769, 0.211358, 0.199668, 2.87505, -4.21223, 3.06169;
% 0.350341, 0.214214, 0.22534, 2.95815, -4.21223, 3.01951;
% 0.416544, 0.215846, 0.256717, 3.03229, -4.21351, 2.98116;
% 0.467533, 0.218294, 0.284833, 3.10005, -4.22245, 2.94664;
% 0.502485, 0.21911, 0.301948, 3.14862, -4.23524, 2.9198;
% 0.516466, 0.217886, 0.304393, 3.17803, -4.24035, 2.89551;
% 0.51441, 0.217886, 0.303985, 3.19081, -4.24163, 2.87633;
% 0.513176, 0.220334, 0.303985, 3.19081, -4.24035, 2.87505;
% 0.513587, 0.240736, 0.303985, 3.19209, -4.24035, 2.8725;
% 0.515643, 0.285211, 0.307652, 3.19209, -4.24163, 2.86483;
% 0.5214, 0.353759, 0.324359, 3.19081, -4.24674, 2.8546;
% 0.527979, 0.441077, 0.347994, 3.19081, -4.26208, 2.83414;
% 0.547717, 0.530027, 0.382222, 3.19848, -4.29149, 2.80474;
% 0.571155, 0.625505, 0.424601, 3.22277, -4.31833, 2.77406;
% 0.59007, 0.727919, 0.475537, 3.25601, -4.34134, 2.72037;
% 0.594182, 0.83319, 0.529325, 3.2803, -4.38353, 2.66412;
% 0.594182, 0.938053, 0.592078, 3.29819, -4.42827, 2.61299;
% 0.597472, 1.0372, 0.661758, 3.3212, -4.47046, 2.55163;
% 0.600762, 1.08821, 0.696802, 3.35061, -4.51264, 2.48515;
% 0.601995, 1.08821, 0.696802, 3.35955, -4.52543, 2.45703;
% 0.60364, 1.08943, 0.696802, 3.36083, -4.52415, 2.45575;
% 0.60364, 1.09596, 0.696802, 3.36083, -4.51137, 2.45575;
% 0.603229, 1.10861, 0.696802, 3.36211, -4.47429, 2.45575;
% 0.59336, 1.13676, 0.706989, 3.36339, -4.43211, 2.44936;
% 0.57979, 1.17471, 0.721659, 3.36339, -4.39887, 2.4289;
% 0.567866, 1.2147, 0.737958, 3.36467, -4.35285, 2.40973;
% 0.55553, 1.25795, 0.748553, 3.36339, -4.29916, 2.39566;
% 0.540315, 1.2914, 0.762408, 3.36211, -4.23396, 2.38416;
% 0.515232, 1.32935, 0.775447, 3.36339, -4.18155, 2.37393;
% 0.508242, 1.34322, 0.781967, 3.36339, -4.14447, 2.37138;
% 0.508653, 1.35057, 0.785634, 3.36339, -4.1253, 2.3701;
% 0.507831, 1.35179, 0.785634, 3.36595, -4.12018, 2.3701;
% 0.501663, 1.35179, 0.785634, 3.3685, -4.11763, 2.37138;
% 0.483159, 1.35261, 0.785634, 3.37106, -4.11763, 2.3701;
% 0.432581, 1.36485, 0.783189, 3.3685, -4.11763, 2.3701;
% 0.351164, 1.37056, 0.773817, 3.357, -4.13041, 2.37265;
% 0.220403, 1.38117, 0.761593, 3.30331, -4.16493, 2.411;
% 0.0476991, 1.38729, 0.757518, 3.21893, -4.20839, 2.49921;
% -0.146387, 1.38729, 0.757518, 3.10644, -4.26592, 2.62449;
% -0.290306, 1.37179, 0.757518, 2.98116, -4.31322, 2.75105;
% -0.427236, 1.3575, 0.757518, 2.89039, -4.32472, 2.85076;
% -0.533736, 1.34771, 0.758333, 2.82647, -4.32728, 2.93514;
% -0.648461, 1.33955, 0.757925, 2.75105, -4.32728, 3.01312;
% -0.703561, 1.32813, 0.758333, 2.66029, -4.32856, 3.08982;
% -0.706851, 1.32527, 0.758333, 2.61682, -4.33239, 3.12178;
% -0.706851, 1.32527, 0.758333, 2.61682, -4.33239, 3.12561;
% -0.702328, 1.31425, 0.753443, 2.61682, -4.33239, 3.12561;
% -0.690403, 1.26121, 0.709434, 2.61682, -4.32984, 3.12689;
% -0.668198, 1.169, 0.627122, 2.6181, -4.30043, 3.12945;
% -0.638181, 1.04822, 0.515878, 2.62833, -4.22629, 3.13968;
% -0.615154, 0.897659, 0.404227, 2.66668, -4.1368, 3.15246;
% -0.600762, 0.727511, 0.314172, 2.69736, -4.04093, 3.1793;
% -0.596238, 0.576134, 0.265274, 2.72037, -3.95911, 3.20359;
% -0.586781, 0.430468, 0.23349, 2.74977, -3.91181, 3.23172;
% -0.584314, 0.299492, 0.21719, 2.77278, -3.88496, 3.25601;
% -0.574034, 0.180348, 0.205781, 2.79452, -3.87729, 3.27774;
% -0.567866, 0.0714047, 0.207003, 2.82647, -3.87602, 3.30714;
% -0.559642, 0.0142809, 0.211078, 2.83798, -3.87602, 3.31993;
% -0.556352, 0.014689, 0.207411, 2.84054, -3.88241, 3.33015;
% -0.550595, 0.015097, 0.204558, 2.84054, -3.88496, 3.33654;
% -0.547717, 0.015097, 0.204558, 2.84309, -3.88496, 3.33654;
% -0.540727, 0.0179532, 0.203336, 2.84309, -3.88496, 3.33654;
% -0.513999, 0.030602, 0.191926, 2.84565, -3.88496, 3.33654;
% -0.466711, 0.0526354, 0.165847, 2.86866, -3.88113, 3.33654;
% -0.394751, 0.0738528, 0.13773, 2.91213, -3.8709, 3.33654;
% -0.310455, 0.0938461, 0.114504, 2.98244, -3.84917, 3.33654;
% -0.203543, 0.109759, 0.101464, 3.06169, -3.83255, 3.33654;
% -0.0900526, 0.127712, 0.0945369, 3.15502, -3.82999, 3.33527;
% 0.0308399, 0.137913, 0.0941294, 3.24194, -3.82872, 3.32887;
% 0.165713, 0.144441, 0.102687, 3.34421, -3.82872, 3.31225;
% 0.309633, 0.161171, 0.127951, 3.45671, -3.82872, 3.27518;
% 0.441628, 0.179124, 0.170737, 3.58071, -3.82872, 3.2266;
% 0.550595, 0.184428, 0.21556, 3.70599, -3.82616, 3.17675;
% 0.593771, 0.186468, 0.230637, 3.79164, -3.81593, 3.12817;
% 0.584725, 0.188916, 0.230637, 3.82616, -3.80698, 3.10899;
% 0.570333, 0.194629, 0.230637, 3.82999, -3.80826, 3.10644;
% 0.564576, 0.210542, 0.230637, 3.82872, -3.80826, 3.10644;
% 0.56581, 0.250528, 0.230637, 3.82616, -3.80826, 3.10516;
% 0.569922, 0.31418, 0.232675, 3.8236, -3.8121, 3.1026;
% 0.571155, 0.40313, 0.242047, 3.82616, -3.81465, 3.08854;
% 0.575678, 0.507585, 0.261199, 3.82744, -3.81338, 3.05914;
% 0.585547, 0.614488, 0.289723, 3.82872, -3.81465, 3.01567;
% 0.598294, 0.735672, 0.330064, 3.82744, -3.81082, 2.97093;
% 0.611864, 0.880929, 0.398522, 3.81977, -3.80571, 2.90957;
% 0.633246, 1.02007, 0.486947, 3.80059, -3.80571, 2.8188;
% 0.651339, 1.16288, 0.599413, 3.77758, -3.80698, 2.71014;
% 0.659152, 1.31425, 0.735513, 3.74946, -3.81977, 2.56824;
% 0.666553, 1.42238, 0.863057, 3.72389, -3.8594, 2.42251;
% 0.666142, 1.45176, 0.889951, 3.70599, -3.88624, 2.31768;
% 0.66532, 1.45257, 0.889543, 3.69577, -3.89391, 2.29084;
% 0.655862, 1.45257, 0.888728, 3.6996, -3.89391, 2.29084;
% 0.619677, 1.45257, 0.888728, 3.69065, -3.90797, 2.29084;
% 0.527568, 1.44686, 0.887098, 3.64335, -3.95783, 2.28956;
% 0.371313, 1.39586, 0.858574, 3.55131, -4.05882, 2.31641;
% 0.218758, 1.28406, 0.816603, 3.40941, -4.17004, 2.41995;
% 0.113491, 1.12819, 0.771372, 3.25217, -4.22501, 2.55674;
% 0.0394751, 0.971511, 0.739996, 3.12817, -4.23396, 2.67179;
% 0.00452319, 0.810749, 0.714732, 3.03868, -4.22757, 2.74722;
% -0.00863518, 0.654475, 0.68417, 2.96965, -4.18922, 2.78301;
% -0.00781278, 0.518194, 0.646274, 2.90957, -4.11507, 2.79835;
% -0.00616798, 0.43414, 0.62264, 2.8725, -4.01919, 2.80219;
% 0.000822398, 0.405578, 0.611637, 2.84437, -3.96039, 2.80091;
% 0.0016448, 0.39293, 0.606748, 2.81241, -3.92715, 2.80219;
% 0.0016448, 0.385177, 0.607155, 2.79707, -3.91948, 2.80219;
% -0.00205599, 0.378241, 0.607563, 2.79963, -3.92332, 2.80219;
% -0.00370079, 0.368448, 0.610008, 2.80219, -3.93738, 2.80219;
% -0.00452319, 0.353759, 0.613675, 2.80219, -3.97062, 2.80219;
% -0.00534559, 0.335398, 0.616935, 2.80219, -4.01664, 2.80091;
% -0.00575678, 0.317037, 0.621825, 2.79963, -4.07672, 2.80091;
% -0.00493439, 0.296635, 0.627122, 2.7894, -4.13297, 2.80091;
% -0.000411199, 0.276234, 0.629567, 2.77406, -4.18538, 2.80219;
% 0.00452319, 0.258281, 0.631197, 2.74722, -4.23779, 2.80219;
% 0.014392, 0.244408, 0.633234, 2.7127, -4.28893, 2.80091;
% 0.0242607, 0.235839, 0.634049, 2.66412, -4.326, 2.80091;
% 0.0345407, 0.229311, 0.634049, 2.60148, -4.34518, 2.80219;
% 0.0448207, 0.224823, 0.634049, 2.54651, -4.36819, 2.80346;
% 0.055923, 0.224007, 0.634049, 2.49793, -4.38225, 2.80219;
% 0.0657918, 0.224415, 0.633642, 2.44936, -4.38609, 2.79068;
% 0.0740158, 0.228495, 0.631197, 2.411, -4.38481, 2.75872;
% 0.0822398, 0.234207, 0.628752, 2.37393, -4.37842, 2.71526;
% 0.0900526, 0.237471, 0.627122, 2.33814, -4.35924, 2.67563;
% 0.0974541, 0.245224, 0.62264, 2.3049, -4.33751, 2.636;
% 0.104033, 0.253793, 0.617342, 2.2755, -4.31833, 2.60915;
% 0.113902, 0.272154, 0.608785, 2.23842, -4.26592, 2.5593;
% 0.119659, 0.281946, 0.602673, 2.20263, -4.22118, 2.50944;
% 0.126238, 0.291739, 0.595746, 2.17451, -4.18922, 2.47109;
% 0.132817, 0.30194, 0.588411, 2.14766, -4.15598, 2.43529;
% 0.136107, 0.316221, 0.581076, 2.12465, -4.11763, 2.39439;
% 0.13734, 0.336622, 0.571704, 2.10803, -4.0601, 2.36115;
% 0.138985, 0.362736, 0.561517, 2.09013, -3.97701, 2.32152;
% 0.138985, 0.384769, 0.549292, 2.07607, -3.89647, 2.28572;
% 0.136518, 0.407619, 0.538697, 2.06584, -3.82488, 2.26144;
% 0.130761, 0.428836, 0.528103, 2.06073, -3.729, 2.23715;
% 0.12336, 0.448421, 0.517916, 2.06201, -3.64591, 2.22564;
% 0.112257, 0.466374, 0.509766, 2.07096, -3.58327, 2.22436;
% 0.0982765, 0.480247, 0.504061, 2.10547, -3.51807, 2.22564;
% 0.0842958, 0.490856, 0.501209, 2.16044, -3.45288, 2.22564;
% 0.066203, 0.498609, 0.499986, 2.23331, -3.40813, 2.23715;
% 0.0505775, 0.499017, 0.499579, 2.3228, -3.36595, 2.26783;
% 0.0365967, 0.499833, 0.499579, 2.40333, -3.34933, 2.29979;
% 0.0172704, 0.499833, 0.499986, 2.49026, -3.34933, 2.34325;
% -0.0012336, 0.499017, 0.499986, 2.60787, -3.34933, 2.3995;
% -0.0193263, 0.496568, 0.499986, 2.73699, -3.34933, 2.4647;
% -0.0378303, 0.493304, 0.499579, 2.84693, -3.34933, 2.56952;
% -0.0551006, 0.486776, 0.494689, 2.93641, -3.35572, 2.68202;
% -0.0740158, 0.474943, 0.488984, 3.03229, -3.39918, 2.7779;
% -0.0970429, 0.452502, 0.486132, 3.12178, -3.46055, 2.86994;
% -0.115547, 0.425572, 0.486132, 3.19081, -3.54492, 2.94792;
% -0.128294, 0.395378, 0.489391, 3.23427, -3.6344, 2.99011;
% -0.136518, 0.363144, 0.493874, 3.25856, -3.73923, 3.00545;
% -0.14063, 0.331726, 0.503246, 3.27007, -3.85045, 3.00672;
% -0.141041, 0.303572, 0.514248, 3.27007, -3.93866, 3.00672;
% -0.138985, 0.28113, 0.526065, 3.27135, -4.03453, 3.00672;
% -0.134462, 0.262769, 0.53829, 3.26623, -4.11124, 3.00672;
% -0.128294, 0.244408, 0.550107, 3.24322, -4.17004, 3.008;
% -0.121304, 0.225639, 0.562332, 3.21382, -4.23012, 3.008;
% -0.111024, 0.206461, 0.573741, 3.18314, -4.29276, 3.00417;
% -0.100333, 0.190548, 0.585151, 3.15118, -4.35668, 3.00033;
% -0.088819, 0.179124, 0.597375, 3.10516, -4.41165, 2.99778;
% -0.074427, 0.174635, 0.605118, 3.04635, -4.45256, 2.98627;
% -0.0608574, 0.176676, 0.607563, 2.98627, -4.47941, 2.9671;
% -0.0472879, 0.18402, 0.60797, 2.91724, -4.48708, 2.95048;
% -0.0378303, 0.196669, 0.607155, 2.8661, -4.48708, 2.93258;
% -0.0283727, 0.21503, 0.59819, 2.80858, -4.45767, 2.89295;
% -0.0230271, 0.242776, 0.585558, 2.74082, -4.38481, 2.84181;
% -0.0230271, 0.2999, 0.548885, 2.67691, -4.29021, 2.80219;
% -0.0275503, 0.378241, 0.484094, 2.62321, -4.19689, 2.78301;
% -0.0370079, 0.449237, 0.400559, 2.59125, -4.10229, 2.76895;
% -0.0407087, 0.499425, 0.302763, 2.56313, -3.99107, 2.74722;
% -0.0407087, 0.529619, 0.199261, 2.55163, -3.88624, 2.7306;
% -0.0407087, 0.535739, 0.0896471, 2.55035, -3.78142, 2.71014;
% -0.0353631, 0.535739, -0.00203743, 2.54779, -3.67276, 2.6948;
% -0.0345407, 0.520234, -0.0696802, 2.54907, -3.61012, 2.69225;
% -0.0328959, 0.469639, -0.131211, 2.54779, -3.56154, 2.69225;
% -0.0324847, 0.421083, -0.179702, 2.54779, -3.50785, 2.68841;
% -0.0316623, 0.371304, -0.230637, 2.55163, -3.45415, 2.67946;
% -0.0209711, 0.33295, -0.276683, 2.56569, -3.40174, 2.66924;
% -0.0201487, 0.317853, -0.311727, 2.59509, -3.35828, 2.66029;
% -0.00411199, 0.290515, -0.342289, 2.60276, -3.34805, 2.63855;
% 0.00287839, 0.268073, -0.373665, 2.60276, -3.3391, 2.60148;
% -0.000411199, 0.267257, -0.37448, 2.60148, -3.33399, 2.57975;
% -0.0016448, 0.267257, -0.37448, 2.60276, -3.33399, 2.57847;
% -0.0016448, 0.267257, -0.37448, 2.6066, -3.33399, 2.57719];

% q=[-0.880377, 1.69861, 0.895656, 2.52861, -3.88241, 2.87122;
% -0.880377, 1.69861, 0.895656, 2.53884, -3.85684, 2.86483;
% -0.880377, 1.69861, 0.895656, 2.53628, -3.83894, 2.86227;
% -0.880377, 1.69861, 0.895656, 2.50305, -3.82872, 2.84309;
% -0.880377, 1.69861, 0.895656, 2.45191, -3.81465, 2.81497;
% -0.880377, 1.69861, 0.895656, 2.411, -3.80187, 2.79324;
% -0.880377, 1.69861, 0.895656, 2.36882, -3.80059, 2.76895;
% -0.880377, 1.69861, 0.895656, 2.31768, -3.80059, 2.72932;
% -0.880377, 1.69861, 0.895656, 2.25504, -3.80826, 2.69225;
% -0.880377, 1.69861, 0.895656, 2.19113, -3.83766, 2.65901;
% -0.880377, 1.69861, 0.895656, 2.14383, -3.87346, 2.62961;
% -0.880377, 1.69861, 0.895656, 2.09397, -3.91309, 2.61682;
% -0.880377, 1.69861, 0.895656, 2.04795, -3.96422, 2.61043;
% -0.880377, 1.69861, 0.895656, 2.02749, -4.02431, 2.60915;
% -0.880377, 1.69861, 0.895656, 2.0096, -4.09462, 2.60915;
% -0.880377, 1.69861, 0.895656, 1.98786, -4.15726, 2.60915;
% -0.880377, 1.69861, 0.895656, 1.97764, -4.21223, 2.60787;
% -0.880377, 1.69861, 0.895656, 1.9738, -4.27103, 2.60532;
% -0.880377, 1.69861, 0.895656, 1.97252, -4.32089, 2.60404;
% -0.880377, 1.69861, 0.895656, 1.97252, -4.36307, 2.60532;
% -0.880377, 1.69861, 0.895656, 1.97252, -4.40654, 2.60532;
% -0.880377, 1.69861, 0.895656, 1.97636, -4.4385, 2.6066;
% -0.880377, 1.69861, 0.895656, 1.97636, -4.48068, 2.60787;
% -0.880377, 1.69861, 0.895656, 1.97636, -4.52543, 2.61427;
% -0.880377, 1.69861, 0.895656, 1.98531, -4.57017, 2.61938;
% -0.880377, 1.69861, 0.895656, 2.00832, -4.61236, 2.62577;
% -0.880377, 1.69861, 0.895656, 2.03261, -4.64559, 2.636;
% -0.880377, 1.69861, 0.895656, 2.06201, -4.68267, 2.64239;
% -0.880377, 1.69861, 0.895656, 2.1042, -4.71079, 2.64367;
% -0.880377, 1.69861, 0.895656, 2.16172, -4.74275, 2.64367;
% -0.880377, 1.69861, 0.895656, 2.2397, -4.7632, 2.64367;
% -0.880377, 1.69861, 0.895656, 2.32791, -4.77215, 2.64239;
% -0.880377, 1.69861, 0.895656, 2.40461, -4.77599, 2.64239;
% -0.880377, 1.69861, 0.895656, 2.48131, -4.77854, 2.64111;
% -0.880377, 1.69861, 0.895656, 2.57208, -4.77982, 2.62961;
% -0.880377, 1.69861, 0.895656, 2.67946, -4.77982, 2.60148;
% -0.880377, 1.69861, 0.895656, 2.78173, -4.76704, 2.56824;
% -0.880377, 1.69861, 0.895656, 2.85971, -4.72869, 2.55035;
% -0.880377, 1.69861, 0.895656, 2.91852, -4.67883, 2.55163;
% -0.880377, 1.69861, 0.895656, 2.98116, -4.63153, 2.55035;
% -0.880377, 1.69861, 0.895656, 3.02079, -4.58168, 2.54907;
% -0.880377, 1.69861, 0.895656, 3.03868, -4.51009, 2.54907;
% -0.880377, 1.69861, 0.895656, 3.05402, -4.43338, 2.55035;
% -0.880377, 1.69861, 0.895656, 3.06169, -4.35668, 2.5529;
% -0.880377, 1.69861, 0.895656, 3.06169, -4.28254, 2.56441;
% -0.880377, 1.69861, 0.895656, 3.04891, -4.21862, 2.58998;
% -0.880377, 1.69861, 0.895656, 3.0259, -4.15981, 2.61299;
% -0.880377, 1.69861, 0.895656, 3.00545, -4.10357, 2.62449;
% -0.880377, 1.69861, 0.895656, 2.97093, -4.04987, 2.6475;
% -0.880377, 1.69861, 0.895656, 2.92874, -3.99618, 2.66796;
% -0.880377, 1.69861, 0.895656, 2.884, -3.9655, 2.67818;
% -0.880377, 1.69861, 0.895656, 2.84309, -3.93354, 2.68458;
% -0.880377, 1.69861, 0.895656, 2.77918, -3.90286, 2.69225;
% -0.880377, 1.69861, 0.895656, 2.6948, -3.88369, 2.69352;
% -0.880377, 1.69861, 0.895656, 2.61043, -3.87729, 2.69352;
% -0.880377, 1.69861, 0.895656, 2.54523, -3.87602, 2.69352;
% -0.880377, 1.69861, 0.895656, 2.51327, -3.87218, 2.69352;
% -0.880377, 1.69861, 0.895656, 2.50816, -3.86962, 2.69352];

q=[0.591715, 0.653659, 0.597783, 3.17419, -4.57401, 2.27038;
0.591715, 0.653659, 0.597783, 3.17419, -4.57401, 2.27038;
0.591715, 0.653659, 0.597783, 3.17419, -4.57401, 2.26911;
0.591715, 0.653659, 0.597783, 3.17419, -4.57528, 2.26783;
0.591715, 0.653659, 0.597783, 3.17419, -4.57401, 2.26655;
0.590893, 0.653659, 0.597783, 3.17419, -4.57401, 2.26527;
0.58637, 0.655291, 0.597375, 3.17419, -4.57401, 2.26655;
0.581846, 0.656923, 0.595746, 3.17419, -4.57401, 2.26783;
0.575678, 0.658963, 0.592893, 3.17419, -4.57401, 2.26783;
0.567454, 0.663451, 0.588818, 3.17036, -4.57528, 2.26655;
0.560875, 0.666307, 0.585151, 3.16141, -4.57401, 2.26655;
0.54895, 0.671204, 0.578224, 3.15885, -4.57401, 2.26655;
0.537437, 0.674468, 0.573741, 3.14735, -4.57401, 2.26783;
0.523867, 0.678548, 0.568444, 3.13584, -4.57401, 2.26783;
0.51112, 0.684261, 0.560294, 3.12561, -4.57401, 2.26783;
0.501251, 0.686709, 0.557442, 3.11155, -4.57401, 2.26783;
0.488504, 0.690789, 0.550922, 3.09877, -4.57401, 2.26655;
0.473701, 0.696501, 0.543995, 3.08215, -4.57401, 2.26527;
0.460132, 0.700174, 0.53829, 3.06809, -4.57401, 2.26783;
0.448207, 0.702214, 0.535845, 3.05786, -4.57401, 2.26783;
0.434226, 0.705886, 0.529733, 3.03613, -4.57401, 2.26911;
0.423946, 0.709558, 0.52525, 3.01951, -4.57273, 2.27422;
0.423535, 0.709558, 0.52525, 3.01312, -4.57273, 2.28828;
0.422712, 0.709558, 0.52525, 3.01312, -4.57145, 2.29339;
0.413666, 0.713231, 0.519138, 3.01056, -4.56889, 2.29212;
0.397629, 0.718535, 0.511396, 2.99522, -4.56506, 2.29212;
0.379537, 0.721799, 0.506506, 2.98755, -4.56506, 2.29723;
0.36021, 0.726695, 0.499579, 2.97093, -4.56506, 2.31257;
0.343762, 0.730776, 0.493059, 2.9492, -4.56122, 2.33047;
0.328548, 0.73404, 0.488169, 2.92874, -4.55867, 2.34581;
0.310455, 0.737304, 0.484502, 2.91213, -4.55867, 2.36243;
0.288662, 0.740568, 0.478797, 2.89679, -4.55739, 2.37265;
0.269335, 0.743424, 0.473907, 2.884, -4.55611, 2.37777;
0.251243, 0.745873, 0.47024, 2.86738, -4.55739, 2.39566;
0.231094, 0.748321, 0.466165, 2.85332, -4.55739, 2.40461;
0.210534, 0.749545, 0.46372, 2.83542, -4.55483, 2.41228;
0.187918, 0.752401, 0.45883, 2.81497, -4.55227, 2.42507;
0.169825, 0.754441, 0.455163, 2.79579, -4.54716, 2.44169;
0.155844, 0.755665, 0.453125, 2.78045, -4.54332, 2.45447;
0.135696, 0.756889, 0.451903, 2.75361, -4.54332, 2.45958;
0.110201, 0.757705, 0.44905, 2.72037, -4.54205, 2.46086;
0.0863518, 0.758521, 0.447828, 2.69225, -4.53949, 2.47492;
0.068259, 0.759338, 0.446198, 2.66796, -4.53821, 2.49154;
0.0501663, 0.759746, 0.445383, 2.63855, -4.53949, 2.49666;
0.0308399, 0.759746, 0.445383, 2.61299, -4.53949, 2.50177;
0.00699038, 0.759746, 0.445383, 2.58998, -4.53949, 2.51583;
-0.0139808, 0.759338, 0.445383, 2.5593, -4.53949, 2.52861;
-0.0333071, 0.758521, 0.445383, 2.54012, -4.53821, 2.54396;
-0.0522223, 0.758521, 0.445383, 2.5235, -4.53821, 2.55674;
-0.068259, 0.757297, 0.44579, 2.50049, -4.53821, 2.56824;
-0.0879966, 0.755257, 0.448235, 2.48131, -4.53821, 2.57464;
-0.11678, 0.753217, 0.45068, 2.45575, -4.53821, 2.57719;
-0.142686, 0.750361, 0.453125, 2.43018, -4.53821, 2.59637;
-0.165713, 0.747913, 0.455978, 2.41484, -4.53821, 2.61554;
-0.184628, 0.745873, 0.458015, 2.39566, -4.53821, 2.6181;
-0.196553, 0.745057, 0.459237, 2.38672, -4.53821, 2.62577;
-0.201487, 0.744649, 0.459237, 2.37777, -4.53949, 2.63216;
-0.201487, 0.744649, 0.459645, 2.37138, -4.53949, 2.64111;
-0.204777, 0.743833, 0.460867, 2.37138, -4.53949, 2.64878;
-0.229449, 0.738936, 0.465757, 2.35603, -4.53949, 2.65773;
-0.255766, 0.734856, 0.47187, 2.32791, -4.53949, 2.67307;
-0.273858, 0.732, 0.475537, 2.30618, -4.53949, 2.68458;
-0.29483, 0.727511, 0.480834, 2.287, -4.53949, 2.68841;
-0.319913, 0.721799, 0.488984, 2.25504, -4.53949, 2.69352;
-0.342117, 0.715271, 0.497134, 2.21797, -4.53949, 2.71142;
-0.361855, 0.71119, 0.502431, 2.19368, -4.53949, 2.72804;
-0.375425, 0.707518, 0.507321, 2.17962, -4.54077, 2.74338;
-0.386116, 0.70303, 0.513433, 2.16811, -4.53949, 2.76383;
-0.406676, 0.696093, 0.522805, 2.14638, -4.53821, 2.78045;
-0.426824, 0.690381, 0.530548, 2.12721, -4.53821, 2.78301;
-0.437927, 0.686301, 0.535845, 2.10931, -4.54077, 2.79579;
-0.444095, 0.683853, 0.539512, 2.09908, -4.5446, 2.80219;
-0.448207, 0.681404, 0.542365, 2.09525, -4.54716, 2.80474;
-0.451085, 0.680588, 0.543587, 2.09269, -4.54972, 2.80986;
-0.456842, 0.677324, 0.54807, 2.09269, -4.55355, 2.81369;
-0.469178, 0.672428, 0.553774, 2.08119, -4.55611, 2.81625;
-0.477813, 0.668756, 0.558664, 2.06329, -4.55867, 2.83287;
-0.488915, 0.663451, 0.565999, 2.0569, -4.56122, 2.84309;
-0.505363, 0.656923, 0.574964, 2.04795, -4.56378, 2.84693;
-0.5214, 0.647538, 0.586781, 2.03772, -4.56761, 2.86099;
-0.539493, 0.637746, 0.599413, 2.02494, -4.57912, 2.884;
-0.556763, 0.629585, 0.6096, 2.00576, -4.58423, 2.90446;
-0.56951, 0.624281, 0.61612, 1.98659, -4.58551, 2.90957;
-0.577323, 0.619384, 0.62264, 1.97892, -4.59062, 2.9134;
-0.590482, 0.612448, 0.630789, 1.96741, -4.59318, 2.91852;
-0.599117, 0.606328, 0.638939, 1.95591, -4.59446, 2.9313;
-0.608574, 0.601023, 0.644644, 1.94696, -4.59702, 2.94408;
-0.617621, 0.594903, 0.652794, 1.93673, -4.59829, 2.94792;
-0.6242, 0.590415, 0.658091, 1.92906, -4.60469, 2.95048;
-0.629546, 0.586334, 0.662981, 1.92139, -4.60596, 2.95431;
-0.632424, 0.584702, 0.665426, 1.91116, -4.60596, 2.95815;
-0.633246, 0.583886, 0.665833, 1.90605, -4.60724, 2.9607;
-0.633658, 0.583478, 0.666648, 1.90221, -4.60596, 2.96198;
-0.636947, 0.581438, 0.669093, 1.90349, -4.60724, 2.96326;
-0.637358, 0.580214, 0.670723, 1.90221, -4.60852, 2.96454;
-0.637358, 0.580214, 0.671131, 1.90349, -4.6098, 2.96326;
-0.649694, 0.572461, 0.681318, 1.86514, -4.60724, 2.96965;
-0.649694, 0.572461, 0.681318, 1.86642, -4.60852, 2.96965;
-0.649694, 0.572461, 0.681318, 1.86514, -4.60724, 2.96965;
-0.649694, 0.572461, 0.681318, 1.86514, -4.60724, 2.96965;
-0.651339, 0.574093, 0.681318, 1.86514, -4.60724, 2.96965;
-0.656685, 0.577358, 0.678873, 1.86258, -4.60852, 2.96837;
-0.659563, 0.579806, 0.67602, 1.86003, -4.60852, 2.96837;
-0.662853, 0.581846, 0.67276, 1.85491, -4.6098, 2.96965;
-0.666965, 0.585518, 0.669501, 1.8498, -4.60469, 2.96965;
-0.671899, 0.590006, 0.665426, 1.84597, -4.59829, 2.96837;
-0.676422, 0.594495, 0.659721, 1.84213, -4.58679, 2.96837;
-0.681357, 0.601023, 0.651164, 1.8383, -4.57528, 2.96837;
-0.688347, 0.608368, 0.641791, 1.83446, -4.5625, 2.96837;
-0.695749, 0.615304, 0.633234, 1.82679, -4.54972, 2.96837;
-0.702739, 0.622649, 0.62264, 1.81912, -4.53949, 2.96837;
-0.708084, 0.628769, 0.614897, 1.81656, -4.53182, 2.96837;
-0.715075, 0.63693, 0.604303, 1.81528, -4.52415, 2.96837;
-0.722888, 0.646722, 0.590041, 1.81273, -4.51137, 2.96582;
-0.7307, 0.656515, 0.577816, 1.81273, -4.49986, 2.96454;
-0.739336, 0.664675, 0.566406, 1.81145, -4.48835, 2.96454;
-0.749204, 0.671204, 0.555812, 1.80889, -4.47685, 2.96326;
-0.757428, 0.67814, 0.54644, 1.8025, -4.4679, 2.96326;
-0.768531, 0.684669, 0.537882, 1.79355, -4.46151, 2.96198;
-0.779633, 0.689565, 0.530955, 1.77949, -4.45256, 2.95687;
-0.791558, 0.693645, 0.524028, 1.77182, -4.44106, 2.95303;
-0.801427, 0.698542, 0.517508, 1.76671, -4.43211, 2.95175;
-0.811707, 0.703846, 0.509766, 1.75904, -4.42444, 2.9492;
-0.82322, 0.707926, 0.502839, 1.75009, -4.41421, 2.94792;
-0.833911, 0.712822, 0.495504, 1.74114, -4.39887, 2.94664;
-0.842547, 0.719351, 0.486947, 1.73858, -4.38736, 2.94664;
-0.854471, 0.725471, 0.477167, 1.73219, -4.37969, 2.94536;
-0.86763, 0.731592, 0.468202, 1.72324, -4.3733, 2.94281;
-0.879966, 0.736896, 0.459645, 1.70918, -4.36819, 2.94153;
-0.889012, 0.741384, 0.45231, 1.69767, -4.36052, 2.93769;
-0.89847, 0.745465, 0.44579, 1.68617, -4.35157, 2.93258;
-0.908338, 0.749545, 0.438456, 1.6785, -4.34262, 2.92747;
-0.917796, 0.754441, 0.430713, 1.67339, -4.33112, 2.92491;
-0.926842, 0.758113, 0.424601, 1.66955, -4.32217, 2.92107;
-0.935066, 0.762194, 0.417266, 1.66572, -4.31705, 2.91852;
-0.944113, 0.767498, 0.409524, 1.65932, -4.3145, 2.91724;
-0.954804, 0.77117, 0.401782, 1.6491, -4.30938, 2.91596;
-0.966317, 0.775659, 0.394855, 1.6312, -4.30299, 2.91596;
-0.978242, 0.778923, 0.388742, 1.62353, -4.2966, 2.91468;
-0.988933, 0.781371, 0.383852, 1.61714, -4.28637, 2.91085;
-0.990167, 0.785451, 0.377333, 1.60947, -4.2787, 2.90829;
-0.990578, 0.788715, 0.37122, 1.60947, -4.27359, 2.89934;
-0.990578, 0.790348, 0.368775, 1.61075, -4.27231, 2.88528;
-0.990578, 0.790348, 0.368368, 1.61075, -4.27103, 2.87633;
-0.990578, 0.790348, 0.368368, 1.60947, -4.27231, 2.87505;
-0.990578, 0.790348, 0.368368, 1.60947, -4.27231, 2.87377;
-0.990578, 0.790348, 0.368368, 1.60819, -4.27231, 2.87505;
-0.990578, 0.790348, 0.368368, 1.60819, -4.27231, 2.87377;
-0.990578, 0.790348, 0.368368, 1.60691, -4.27231, 2.87377;
-0.990578, 0.790348, 0.368368, 1.60691, -4.27231, 2.87377;
-0.990578, 0.790348, 0.368368, 1.60819, -4.27231, 2.87377;
-0.990578, 0.790348, 0.368775, 1.61075, -4.27359, 2.87377;
-0.990578, 0.790348, 0.368775, 1.60947, -4.27359, 2.87377;
-0.990578, 0.790348, 0.368775, 1.60947, -4.27487, 2.87377;
-0.990578, 0.790348, 0.367553, 1.60947, -4.27359, 2.87377;
-0.990578, 0.790348, 0.365923, 1.60947, -4.27487, 2.87377;
-0.990578, 0.790348, 0.365923, 1.60819, -4.27487, 2.87377;
-0.990578, 0.789123, 0.365923, 1.6133, -4.27359, 2.87377;
-0.990578, 0.788307, 0.365923, 1.61714, -4.26975, 2.87377;
-0.990578, 0.787083, 0.36796, 1.62097, -4.25569, 2.87377;
-0.990578, 0.79198, 0.371628, 1.63503, -4.24035, 2.87377;
-0.990578, 0.800956, 0.37448, 1.6606, -4.21862, 2.87377;
-0.990578, 0.80422, 0.37774, 1.67594, -4.19944, 2.87377;
-0.990578, 0.809933, 0.38263, 1.70151, -4.18027, 2.87377;
-0.990578, 0.821766, 0.393632, 1.73219, -4.15854, 2.87377;
-0.990578, 0.830742, 0.403819, 1.77054, -4.13552, 2.87377;
-0.990578, 0.836863, 0.411562, 1.81145, -4.12018, 2.87377;
-0.990578, 0.842575, 0.422156, 1.85875, -4.11507, 2.87377;
-0.990989, 0.855632, 0.436826, 1.91628, -4.11124, 2.87377;
-0.990989, 0.86012, 0.447828, 1.95846, -4.10612, 2.87377;
-0.990989, 0.862976, 0.452718, 2.00704, -4.09973, 2.87377;
-0.990989, 0.872769, 0.45883, 2.05562, -4.09078, 2.87377;
-0.990989, 0.874401, 0.463312, 2.10803, -4.08567, 2.87377;
-0.990989, 0.874401, 0.466572, 2.15661, -4.08183, 2.87377;
-0.991401, 0.874401, 0.466572, 2.20263, -4.08183, 2.8725;
-0.991401, 0.874809, 0.467795, 2.26271, -4.08311, 2.87377;
-0.991401, 0.874809, 0.467795, 2.31257, -4.07288, 2.8725;
-0.991401, 0.875217, 0.472685, 2.35731, -4.06905, 2.87377;
-0.991401, 0.876033, 0.476352, 2.40206, -4.05882, 2.87377;
-0.991401, 0.878073, 0.477982, 2.43913, -4.05627, 2.87377;
-0.991812, 0.88297, 0.481649, 2.48131, -4.05882, 2.87377;
-0.991812, 0.88501, 0.483687, 2.53245, -4.08567, 2.87377;
-0.991812, 0.885418, 0.487762, 2.58103, -4.13169, 2.87377;
-0.991812, 0.884602, 0.491836, 2.64239, -4.17515, 2.8725;
-0.991812, 0.884602, 0.491836, 2.69225, -4.20456, 2.86099;
-0.991812, 0.883786, 0.490206, 2.72548, -4.23012, 2.85588;
-0.991401, 0.877665, 0.486539, 2.75744, -4.25186, 2.85588;
-0.991401, 0.863792, 0.486947, 2.80986, -4.28382, 2.86355;
-0.991401, 0.845431, 0.486539, 2.8661, -4.30299, 2.90318;
-0.991401, 0.844615, 0.484094, 2.90701, -4.30427, 2.94281;
-0.991812, 0.844207, 0.474314, 2.92491, -4.30427, 2.93514;
-0.991812, 0.842575, 0.467387, 2.94281, -4.30555, 2.86483;
-0.991812, 0.83523, 0.456793, 2.97732, -4.31961, 2.79835;
-0.991812, 0.829518, 0.449458, 2.99778, -4.35285, 2.7715;
-0.991812, 0.827886, 0.448235, 3.01567, -4.37842, 2.74849;
-0.991812, 0.818909, 0.444975, 3.03613, -4.40143, 2.74722;
-0.991812, 0.818501, 0.444568, 3.05914, -4.42188, 2.74722;
-0.991812, 0.816869, 0.442938, 3.07959, -4.43083, 2.74466;
-0.991812, 0.816461, 0.440086, 3.09365, -4.44745, 2.73955;
-0.991812, 0.814829, 0.438456, 3.11027, -4.46023, 2.73188;
-0.991812, 0.813197, 0.436826, 3.12561, -4.47301, 2.72293;
-0.991812, 0.809525, 0.433566, 3.12945, -4.47685, 2.72293;
-0.991812, 0.806669, 0.430306, 3.12945, -4.47685, 2.72293;
-0.991812, 0.806261, 0.430306, 3.13073, -4.47685, 2.72421;
-0.991812, 0.805853, 0.429898, 3.12945, -4.47685, 2.72293];

for i=1:size(q,1)
% Forward kinematics of the Phantom
[stylus_pose(:,i),T_phantom(:,:,i),T_base(:,:,i),T_link1(:,:,i),T_link2(:,:,i)]=Phantom_FK(q(i,:),d1,L1,L2);
end

% Plot of the Phantom
scrsz = get(0,'ScreenSize');
figure('Position',[0 0 scrsz(3) scrsz(4)])
% figure
axis equal
axis([-300 300 -300 300 0 600]);
% Verteces defining the base
[a_base,z_base]=ndgrid(linspace(0,1,no_arc)*2*pi,linspace(0,d1,no_slice));
y_base=cos(a_base)*r_base;
x_base=sin(a_base)*r_base;
% Verteces defining the links
[a_link,x_link]=ndgrid(linspace(0,1,no_arc)*2*pi,linspace(0,L1,no_slice));
y_link=cos(a_link)*r_link;
z_link=sin(a_link)*r_link;
% Verteces defining the stylus
[a_stylus,x_stylus]=ndgrid(linspace(0,1,no_arc)*2*pi,linspace(-41,135,no_slice));
y_stylus=cos(a_stylus)*r_stylus;
z_stylus=sin(a_stylus)*r_stylus;
r_pt=size(q,1);
F(r_pt)=struct('cdata',[],'colormap',[]);
% Prepare the new file.
vidObj = VideoWriter('sim2.avi','Uncompressed AVI');
open(vidObj);
for i=1:r_pt
    view(3)
    grid on
    hold on
    % Plot of the base
    surf(x_base,y_base,z_base,z_base*0);
    colormap bone
    shading interp
    freezeColors
    % Plot of the first link
    Link1=T_link1(1:3,1:3,i)*[reshape(x_link,1,no_arc*no_slice);reshape(y_link,1,no_arc*no_slice);...
        reshape(z_link,1,no_arc*no_slice)];
    surf(reshape(Link1(1,:),no_arc,no_slice)+T_base(1,4,i),reshape(Link1(2,:),no_arc,no_slice)+T_base(2,4,i),...
        reshape(Link1(3,:),no_arc,no_slice)+T_base(3,4,i),z_link*0);
    colormap spring
    shading interp
    freezeColors
    % Plot of the second link
    Link2=T_link2(1:3,1:3,i)*[reshape(x_link,1,no_arc*no_slice);reshape(y_link,1,no_arc*no_slice);...
        reshape(z_link,1,no_arc*no_slice)];
    surf(reshape(Link2(1,:),no_arc,no_slice)+T_link1(1,4,i),reshape(Link2(2,:),no_arc,no_slice)+T_link1(2,4,i),...
        reshape(Link2(3,:),no_arc,no_slice)+T_link1(3,4,i),z_link*0);
    colormap winter
    shading interp
    freezeColors
    % Plot of the stylus
    T_stylus=T_phantom(:,:,i)*[reshape(x_stylus,1,no_arc*no_slice);reshape(y_stylus,1,no_arc*no_slice);...
        reshape(z_stylus,1,no_arc*no_slice);ones(1,no_arc*no_slice)];
    surf(reshape(T_stylus(1,:),no_arc,no_slice),reshape(T_stylus(2,:),no_arc,no_slice),...
        reshape(T_stylus(3,:),no_arc,no_slice),z_stylus*0);
    colormap gray
    shading interp
    freezeColors
    F(i)=getframe;
    writeVideo(vidObj,F(i));
    if i==r_pt
    else clf
    end
    axis equal
    axis([-300 300 -300 300 0 600]);
end
% Close the file.
close(vidObj);
