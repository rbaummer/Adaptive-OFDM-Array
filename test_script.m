%Author: Rob Baummer
%
function test_script

close all;
%% Network Test
%create first node in network
%128 channels 
%network = adhoc_network(AccessPoint([0 0 0], AntennaConfiguration('Linear', 10, 0, 5e6)));
%16 channels
network = adhoc_network(AccessPoint([0 0 0], AntennaConfiguration('Linear', 10, 0, 600e3)));
%add additional nodes
network.nodes = AccessPoint([5000 10000 0], AntennaConfiguration('Circular', 10, 0, 5e6));
network.nodes = AccessPoint([10000 5000 0]);
network.nodes = AccessPoint([-12000 7000 0]);
% network.nodes = AccessPoint([0 0 10000]);
% network.nodes = AccessPoint([10000 10000 10000]);

%Calculate AoA for node 1
[d,az,el] = network.CalculateEnvironment;
angle = az(1,:)*180/pi %#ok<NOPRT>

%Calculate input at node 1 array
array_waveform = network.CalculateArrayInput(1);

%CMA algorithm
bf = CMABeamformer(network.nodes(1).OFDM_inst);
%128 channelsbf.FrequencyDomainCMA(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.0005);
bf.FrequencyDomainCMA(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.001);
%bf.TimeDomainCMA(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.00001);

%K-Omega plot of input at node 1 array
L = 1024;
M = L/2;
twin = ones(L,1);
%swin = ones(size(array_waveform,1),1);
swin = hamming(size(array_waveform,1));
% twin = hamming(L);
% swin = hamming(nrcv);
ntFFT = 2048;
nsFFT = 2048;
Ts = network.nodes(1).OFDM_inst.Ts/(network.nodes(1).OFDM_inst.L*network.nodes(1).OFDM_inst.N);
network.nodes(1).AntennaConfiguration.PlotKOmega(array_waveform(:,1:50000), 1/Ts, L, M, twin, swin, ntFFT, nsFFT)

%demodulate input of 1 antenna
network.nodes(1).OFDM_inst.demodulate(array_waveform(1,:));
%demodulate input of all antennas 
network.nodes(1).OFDM_inst.demodulate(sum(array_waveform,1));

% [d,az,el] = network.CalculateEnvironment;
% 
% network.PlotEnvironment(1);
% 
% lambda = 3e8/800e6;
% V = network.nodes(1).AntennaConfiguration.CalculateReplicaVectors(az(1,:), el(1,:), lambda);
% network.nodes(3).AntennaConfiguration.PlotBeampattern2D(0);
% network.nodes(1).AntennaConfiguration.PlotBeampattern3D;
% network.nodes(2).AntennaConfiguration.PlotBeampattern3D;

temp = 1;

