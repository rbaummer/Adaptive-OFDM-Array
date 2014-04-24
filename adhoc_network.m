%Author: Rob Baummer
%
%adhoc_network is constructed of nodes.  Nodes consist of AccessPoints
%
%properties
%   nodes: Array of AccessPoints in the network
%
%default constructor: adhoc_network(new_AccessPoint)
%   if new_AccessPoint is provided it becomes the first node
%   if new_AccessPoint is not provided the default AccessPoint([0 0 0])
%     becomes the first node
%method set: obj.nodes = new_AccessPoint
%   if new_AccessPoint is an array of access points the nodes are set to
%     the array
%   if new_AccessPoint is a single access point a new node at the end of
%     the adhoc_network is added for the access point
%method calculate_environment: [d, az, el] = obj.calculate_environment
%   calculates the distance between each node and all other nodes (d)
%   calculates angle of arrival at each node for signals originating at all
%   other nodes (azimuth, elevation)
classdef adhoc_network
    %adhoc_network properties
    properties
        nodes
    end
    
    %hidden property for default AccessPoint
    properties(Constant = true, Hidden = true)
        default_AccessPoint = AccessPoint([0 0 0]);
    end
    
    methods
        %constructor
        function obj = adhoc_network(new_AccessPoint)
            if nargin < 1
                obj.nodes = obj.default_AccessPoint;
            else
                obj.nodes = new_AccessPoint;
            end
        end
        
        %set nodes method
        %if new_AccessPoint is a single AccessPoint it is appended to nodes
        %if new_AccessPoint is an array of AccessPoints it replaces nodes
        function obj = set.nodes(obj, new_AccessPoint)
            if nargin < 2
                obj.nodes = [obj.nodes obj.default_AccessPoint];
            else
                if length(new_AccessPoint) == 1
                    obj.nodes = [obj.nodes new_AccessPoint];
                else
                    obj.nodes = new_AccessPoint;
                end
            end
        end
        
        %find the distance between nodes and angle of arrival at each node
        function [d, az, el] = CalculateEnvironment(obj)
            num = length(obj.nodes);
            
            %preallocate matrices for speed
            d = zeros(num,num);
            az = zeros(num,num);
            el = zeros(num,num);
            %for each node in the network
            for i = 1:num
                %calculate environment from all other nodes
                for j = 1:num
                    if i ~= j
                        %deltas from node i to node j
                        delta_x = obj.nodes(j).location(1) - obj.nodes(i).location(1);
                        delta_y = obj.nodes(j).location(2) - obj.nodes(i).location(2);
                        delta_z = obj.nodes(j).location(3) - obj.nodes(i).location(3);
                        
                        %distance between i and j
                        d(i,j) = sqrt(delta_x^2 + delta_y^2 + delta_z^2);
                        
                        %angle of arrival for signal from node j arriving
                        %at node j
                        az(i,j) = atan2(delta_y, delta_x);
                        el(i,j) = atan2(delta_z, sqrt(delta_x^2 + delta_y^2));
                    end
                end
            end
        end 
        
        %Calculate the array inputs at the desired node
        function array_waveform = CalculateArrayInput(obj, desired_node)
            num = length(obj.nodes) - 1;
            
            %calculate the distance and angles of arrival for all signals
            %from other nodes
            [d, az, el] = obj.CalculateEnvironment();
            
            %list of nodes
            node_list = 1:num+1;
            %remove desired node from list
            node_list(desired_node) = [];
            
            %for each node in the network
            %calculate the input at the antenna array for each node
            for i = 1:num
                %current node
                j = node_list(i);
                
                %Perform the broadband simulation for waveform arriving at
                %desired_node from current node.  Result is returned in
                %frequency domain so multiple signals can be summed
                waveform_fft = obj.nodes(desired_node).AntennaConfiguration.BroadbandSimulation(az(desired_node,j),el(desired_node,j), obj.nodes(j).OFDM_inst);

                if i == 1
                    cumulative_fft = waveform_fft;
                else
                    cumulative_fft = cumulative_fft + waveform_fft;
                end
            end
            
            %transform back to time domain and take the real portion
            n=sqrt(1000)*randn(size(cumulative_fft));
            
            %cumulative_fft is a Nxlength matrix.  Take IFFT across rows
            %array_waveform = real(ifft(cumulative_fft,[],2)) + n;
            array_waveform = (ifft(cumulative_fft,[],2)) + n;
        end   
        
        %PlotEnvironment(K) plots the nodes in the network
        %if K is specified it also plots the antenna array at node K
        function PlotEnvironment(obj, K)
            figure;            
            hold on;
            %plot the locaiton of each node in the network
            for i = 1:length(obj.nodes)
                plot3(obj.nodes(i).location(1), obj.nodes(i).location(2), obj.nodes(i).location(3), '*');
            end
            title('Ad-Hoc Network');   
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
            %if a specific node is identified plots its antenna
            %configuration
            if nargin > 1
                N = obj.nodes(K).AntennaConfiguration.N;
                M = obj.nodes(K).AntennaConfiguration.M;
                if M > 0
                    NumElements = M*N;
                else
                    NumElements = N;
                end
                x = ones(NumElements,1)*obj.nodes(K).location(1) + obj.nodes(K).AntennaConfiguration.G(:,1);
                y = ones(NumElements,1)*obj.nodes(K).location(2) + obj.nodes(K).AntennaConfiguration.G(:,2);
                z = ones(NumElements,1)*obj.nodes(K).location(3) + obj.nodes(K).AntennaConfiguration.G(:,3);
                plot3(x, y, z, '.');
            end
        end

    end
    
end