%Author: Rob Baummer
%
function test_script(network)

close all;
%Moved Network setup to separate script to avoid recreating waveforms each
%time the CMA beamformers are rerun
% %% Network Test
% %create first node in network
% %16 channels at location 0,0,0 with 10 element linear array
% network = adhoc_network(AccessPoint([0 0 0], AntennaConfiguration('Linear', 10, 0, 600e3)));
% %add additional node at location -5000,10000,0 with 10 element circular array
% network.nodes = AccessPoint([-5000 10000 0], AntennaConfiguration('Circular', 10, 0, 5e6));
% %add additional nodes at specified location with default array
% network.nodes = AccessPoint([10000 5000 0]);
% network.nodes = AccessPoint([7000 7000 0]);
% network.nodes = AccessPoint([2000 6000 0]);
% network.nodes = AccessPoint([-12000 7000 0]);

%Calculate AoA for node 1
[d,az,el] = network.CalculateEnvironment;
angle = az(1,:)*180/pi %#ok<NOPRT>

%Calculate the broadband input at node 1 array
array_waveform = network.CalculateArrayInput(1);

%CMA algorithm
bf = CMABeamformer(network.nodes(1).OFDM_inst);
%128 channels bf.CMA_FrequencyDomain(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.0005);
[waveform_cma, ~] = bf.CMA_FrequencyDomain(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.01);
[waveform_ocma, ~] = bf.oCMA_FrequencyDomain(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.01);
[waveform_lccma, ~] = bf.LCCMA_FrequencyDomain(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.001, az(1,5));
[waveform_lcocma, ~] = bf.LCoCMA_FrequencyDomain(network.nodes(1).AntennaConfiguration, network.nodes(1).OFDM_inst, array_waveform, 0.005, az(1,5));

%% 2 bits per symbol
%demodulation prior to beamforming from first antenna input
network.nodes(1).OFDM_inst.demodulate(array_waveform(1,:)/mean(abs(array_waveform(1,:))));
title('OFDM QPSK Constellation Diagram: Prior to Beam Forming');
%demodulate
network.nodes(1).OFDM_inst.demodulate(waveform_cma);
title('OFDM QPSK Constellation Diagram: CMA');
%demodulate
network.nodes(1).OFDM_inst.demodulate(waveform_lccma);
title('OFDM QPSK Constellation Diagram: LCCMA');
%demodulate
network.nodes(1).OFDM_inst.demodulate(waveform_ocma);
title('OFDM QPSK Constellation Diagram: oCMA');
%demodulate
network.nodes(1).OFDM_inst.demodulate(waveform_lcocma);
title('OFDM QPSK Constellation Diagram: LCoCMA');

%% 3 bits per symbol
% %demodulation prior to beamforming from first antenna input
% network.nodes(1).OFDM_inst.demodulate(array_waveform(1,:)/mean(abs(array_waveform(1,:))));
% title('OFDM 8QAM Constellation Diagram: Prior to Beam Forming');
% %demodulate
% network.nodes(1).OFDM_inst.demodulate(waveform_cma);
% title('OFDM 8QAM Constellation Diagram: CMA');
% %demodulate
% network.nodes(1).OFDM_inst.demodulate(waveform_lccma);
% title('OFDM 8QAM Constellation Diagram: LCCMA');
% %demodulate
% network.nodes(1).OFDM_inst.demodulate(waveform_ocma);
% title('OFDM 8QAM Constellation Diagram: oCMA');
% %demodulate
% network.nodes(1).OFDM_inst.demodulate(waveform_lcocma);
% title('OFDM 8QAM Constellation Diagram: LCoCMA');
