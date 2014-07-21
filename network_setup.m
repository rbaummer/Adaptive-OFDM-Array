function network = network_setup

%% Network Test
%create first node in network
%16 channels at location 0,0,0 with 10 element linear array
network = adhoc_network(AccessPoint([0 0 0], AntennaConfiguration('Linear', 10, 0, 600e3)));
%add additional node at location -5000,10000,0 with 10 element circular array
network.nodes = AccessPoint([-5000 10000 0], AntennaConfiguration('Circular', 10, 0, 5e6));
%add additional nodes at specified location with default array
network.nodes = AccessPoint([10000 5000 0]);
network.nodes = AccessPoint([7000 7000 0]);
network.nodes = AccessPoint([2000 6000 0]);
network.nodes = AccessPoint([-12000 7000 0]);