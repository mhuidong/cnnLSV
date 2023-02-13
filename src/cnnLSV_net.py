import torch
import torch.nn as nn
from torchvision.transforms import transforms as tt

class Net(nn.Module):
    def __init__(self, line, dropout):
        # l1_in, l1_out = 900, 128
        # l2_in, l2_out = 128, 1
        # dropout = 0.5
        l1_in, l1_out, l2_in, l2_out = line
        super(Net, self).__init__()
        self.conv = nn.Sequential(nn.Conv2d(3, 3, 3, 1, padding=1), nn.ReLU(), nn.MaxPool2d(2),  # 3*100*50
                                  nn.Conv2d(3, 3, 3, 1, padding=1), nn.ReLU(), nn.MaxPool2d(2),  # 3*50*25
                                  nn.Conv2d(3, 3, 3, 1, padding=1), nn.ReLU(), nn.MaxPool2d(2))  # 3*25*12
        self.linear = nn.Sequential(nn.Linear(l1_in, l1_out), nn.ReLU(), nn.Dropout(dropout), nn.Linear(l2_in, l2_out))
        # self.linear = nn.Sequential(nn.Linear(l1_in, l1_out))
        self.sigmoid = nn.Sigmoid()
    def forward(self, x):
        x = self.conv(x)
        x = x.view(x.size(0), -1)
        x = self.linear(x)
        return self.sigmoid(x)


# CPU训练的模型
def pre_sv(img, net_model_file_path, means, std, line, dropout, svtype):
    trans = tt.Compose([tt.ToTensor(), tt.Normalize(means, std)])
    net = Net(line, dropout)
    net.load_state_dict(torch.load(net_model_file_path, map_location=torch.device('cpu'))[svtype])
    img = trans(img)
    x = torch.unsqueeze(img, 0)
    pre = net(x).item()
    if pre >= 0.5:
        return 1
    else:
        return 0