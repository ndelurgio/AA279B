function LLA = eci2lla_datetime(pos,t)
LLA = zeros(length(t),3);
for i = 1:length(t)
    t_vec = datetime2vec(t(i));
    LLA(i,:) = eci2lla(pos(i,:),t_vec);
end
end