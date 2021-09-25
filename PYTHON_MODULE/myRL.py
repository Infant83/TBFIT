import random
import collections
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
import sys
from torch.utils.data import Dataset, DataLoader

class mytb_dataset(Dataset):
    def __init__(self, filename, device=None) -> None:
        if device is None:
            self.device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
        else:
            self.device = device
        self.filename = filename
        data = torch.load(filename, map_location=device)
        self.bands= data['TBA']
        self.target=data['DFT']
        self.n_samples = self.bands.shape[0]
    def __getitem__(self, index):
        return self.bands[index], self.target[0]
    def __len__(self):
        return self.n_samples

class mydqn:
    def __init__(self, buffer_limit=200000, batch_size=64, N_BANDS=13, N_KPOINTS=60, N_W=3, gamma=0.99):
        self.buffer_limit = buffer_limit
        self.batch_size   = batch_size
        self.N_BANDS      = N_BANDS
        self.N_KPOINTS    = N_KPOINTS
        self.N_W          = N_W
        self.gamma        = gamma
        self.device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
        print('Using {} device'.format(self.device))
        self.QNet = self._QNet(self.N_BANDS,self.N_KPOINTS,self.N_W,self.device)
        self.ReplayBuffer = self._ReplayBuffer(self.buffer_limit,self.device)

    class _ReplayBuffer():
        def __init__(self, buffer_limit, device):
            self.buffer = collections.deque(maxlen=buffer_limit)
            self.device = device
        def put(self, transition): 
            self.buffer.append(transition)
        def sample(self, n):
            mini_batch = random.sample(self.buffer, n)
            s_lst, a_lst, r_lst, s_prime_lst, done_mask_lst = [], [], [], [], []
            for transition in mini_batch:
                s, a, r, s_prime, done = transition
                s_lst.append(s)
                a_lst.append([a])
                r_lst.append([r])
                s_prime_lst.append(s_prime)
                done_mask = 0.0 if done else 1.0
                done_mask_lst.append([done_mask])
            return torch.tensor(s_lst,         dtype=torch.float, device=self.device), \
                   torch.tensor(a_lst                           , device=self.device), \
                   torch.tensor(r_lst                           , device=self.device), \
                   torch.tensor(s_prime_lst,   dtype=torch.float, device=self.device), \
                   torch.tensor(done_mask_lst                   , device=self.device)
        def size(self): return len(self.buffer)

    class _QNet(nn.Module):
        def __init__(self,N_BANDS,N_KPOINTS,N_W,device):
            self.N_BANDS = N_BANDS
            self.N_KPOINTS = N_KPOINTS
            self.N_W = N_W
            
            super(mydqn._QNet, self).__init__()
            self.fc1 = nn.Linear(N_KPOINTS*N_BANDS, N_KPOINTS*N_BANDS).to(device)
            self.fc2 = nn.Linear(N_KPOINTS*N_BANDS, N_KPOINTS*N_BANDS).to(device)
            self.fc3 = nn.Linear(N_KPOINTS*N_BANDS, N_KPOINTS*N_BANDS*N_W).to(device)
 
        def forward(self, x):
            x = F.relu(self.fc1(x))
            x = F.relu(self.fc2(x))
            x = self.fc3(x)
            return x
 
        def sample_action(self, obs, epsilon):
            out = self.forward(obs)
            coin = random.random()
            if coin < epsilon:
                return random.randint(0,self.N_KPOINTS*self.N_BANDS*self.N_W-1)
            else :
                return out.argmax().item()


    def train(self, q, q_target, memory, optimizer):
        for i in range(5):
            s,a,r,s_prime,done_mask = memory.sample(self.batch_size)
            q_out = q(s)
            q_a = q_out.gather(1,a)
            max_q_prime = q_target(s_prime).max(1)[0].unsqueeze(1)
            target = r + self.gamma * max_q_prime * done_mask
            loss = F.smooth_l1_loss(q_a, target)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

    def Reporter(self, string=None, step=None, every=None, filenm=None):
        if string is not None:
            if every and step is not None:
                if(step%every == 0):
                    print(string)
            else:
                print(string)
        if filenm is not None:
            print(string, file=open(filenm,"a"))
        else:
            pass
        sys.stdout.flush()

class myddpg:
    def __init__(self, buffer_limit=200000, batch_size=64, N_BANDS=13, N_KPOINTS=60, N_W=3, gamma=0.99, tau=0.005):
        self.buffer_limit = buffer_limit
        self.batch_size   = batch_size
        self.N_BANDS      = N_BANDS
        self.N_KPOINTS    = N_KPOINTS
        self.N_W          = N_W
        self.gamma        = gamma
        self.tau          = tau
        self.device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
        print('Using {} device'.format(self.device))
        self.MuNet = self._MuNet(self.N_BANDS,self.N_KPOINTS,self.N_W,self.device)
        self.QNet = self._QNet(self.N_BANDS,self.N_KPOINTS,self.N_W,self.device)
        self.ReplayBuffer = self._ReplayBuffer(self.buffer_limit,self.device)

    class _ReplayBuffer():
        def __init__(self, buffer_limit, device):
            self.buffer = collections.deque(maxlen=buffer_limit)
            self.device = device

        def put(self, transition):
            self.buffer.append(transition)
 
        def sample(self, n):
            mini_batch = random.sample(self.buffer, n)
            s_lst, a_lst, r_lst, s_prime_lst, done_mask_lst = [], [], [], [], []
 
            for transition in mini_batch:
                s, a, r, s_prime, done = transition
                s_lst.append(s)
                a_lst.append([a])
                r_lst.append([r])
                s_prime_lst.append(s_prime)
                done_mask = 0.0 if done else 1.0
                done_mask_lst.append([done_mask])
 
            return torch.tensor(s_lst,         dtype=torch.float, device=self.device), \
                   torch.tensor(a_lst,         dtype=torch.float, device=self.device), \
                   torch.tensor(r_lst,         dtype=torch.float, device=self.device), \
                   torch.tensor(s_prime_lst,   dtype=torch.float, device=self.device), \
                   torch.tensor(done_mask_lst, dtype=torch.float, device=self.device)
 
        def size(self):
            return len(self.buffer)

    class _MuNet(nn.Module):
        def __init__(self,N_BANDS,N_KPOINTS,N_W,device):
            super(myddpg._MuNet, self).__init__()
            self.N_BANDS=N_BANDS
            self.N_KPOINTS=N_KPOINTS
            self.N_W=N_W

            self.fc1 = nn.Linear(N_BANDS*N_KPOINTS, int(N_BANDS*N_KPOINTS*0.5)).to(device)
            self.fc2 = nn.Linear(int(N_BANDS*N_KPOINTS*0.5), int(N_BANDS*self.N_KPOINTS*0.25)).to(device)
            self.fc_mu = nn.Linear(int(N_BANDS*N_KPOINTS*0.25), 1).to(device)

        def forward(self, x):
            x = F.relu(self.fc1(x))
            x = F.relu(self.fc2(x))
            mu = torch.sigmoid(self.fc_mu(x))*self.N_BANDS*self.N_KPOINTS*self.N_W
            return mu


    class _QNet(nn.Module):
        def __init__(self,N_BANDS,N_KPOINTS,N_W,device):
            super(myddpg._QNet, self).__init__()
            self.fc_s = nn.Linear(N_BANDS*N_KPOINTS, N_W*N_KPOINTS).to(device)
            self.fc_a = nn.Linear(1, N_W*N_KPOINTS).to(device)
            self.fc_q = nn.Linear(N_W*N_KPOINTS*2,  128).to(device)
            self.fc_out = nn.Linear(128,1).to(device)

        def forward(self, x, a):
            h1 = F.relu(self.fc_s(x))
            h2 = F.relu(self.fc_a(a))
            cat = torch.cat([h1,h2], dim=1)
            q = F.relu(self.fc_q(cat))
            q = self.fc_out(q)
            return q

    class OrnsteinUhlenbeckNoise:
        def __init__(self, mu):
            self.theta, self.dt, self.sigma = 0.4, 0.01, 0.2
            self.mu = mu
            self.x_prev = np.zeros_like(self.mu)

        def __call__(self):
            x = self.x_prev + self.theta * (self.mu - self.x_prev) * self.dt + \
                    self.sigma * np.sqrt(self.dt) * np.random.normal(size=self.mu.shape)
            self.x_prev = x
            return x

    def train(self, mu, mu_target, q, q_target, memory, q_optimizer, mu_optimizer):
        s,a,r,s_prime,done_mask  = memory.sample(self.batch_size)

        target = r + self.gamma * q_target(s_prime, mu_target(s_prime)) * done_mask
        q_loss = F.smooth_l1_loss(q(s,a), target.detach())
        q_optimizer.zero_grad()
        q_loss.backward()
        q_optimizer.step()

        mu_loss = -q(s,mu(s)).mean() # That's all for the policy loss.
        mu_optimizer.zero_grad()
        mu_loss.backward()
        mu_optimizer.step()

    def soft_update(self, net, net_target):
        for param_target, param in zip(net_target.parameters(), net.parameters()):
            param_target.data.copy_(param_target.data * (1.0 - self.tau) + param.data * self.tau)

    def Reporter(self, string=None, step=None, every=None, filenm=None):
        if string is not None:
            if every and step is not None:
                if(step%every == 0):
                    print(string)
            else:
                print(string)
        if filenm is not None:
            print(string, file=open(filenm,"a"))
        else:
            pass
        sys.stdout.flush()

def nout_cnn(in_features, network, pooling=None):
    if type(in_features) is int:
        padding = network.padding[0]
        kernel_size = network.kernel_size[0]
        stride = network.stride[0]
        if pooling is not None:
            pool = pooling.kernel_size[0] if type(pooling.kernel_size) is tuple else pooling.kernel_size
        else:
            pool = 1
        n = int( ((in_features + 2 * padding - kernel_size)/stride + 1)/pool)
        n_= float((in_features + 2 * padding - kernel_size)/stride + 1)/pool 
        if abs(n_ - n) > 0.001:
            print("nout is not integer...")
            quit()
        else:
            return n
    elif type(in_features) is list and len(in_features) == 2:
        n = [0,0]
        for i in range(2):
            padding = network.padding[i]
            kernel_size = network.kernel_size[i]
            stride = network.stride[i]
            if pooling is not None:
                pool = pooling.kernel_size[i] if type(pooling.kernel_size) is tuple else pooling.kernel_size
            else:
                pool = 1
            n[i] = int(((in_features[i] + 2 * padding - kernel_size)/stride+1)/pool)
            n_   =     ((in_features[i] + 2 * padding - kernel_size)/stride+1)/pool 
            if abs(n_ - n[i]) > 0.001:
                print("nout is not integer...")
                quit()
    else:
        print("check in_featrues, kernel_size, stride, padding, pooling, etc")
        quit()
    return int(n[0]), int(n[1])

def count_parameters(model):
    nparams = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(model)
    print("nparams: {nparams}")
    return nparams

def get_device():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    return device