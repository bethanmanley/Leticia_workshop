FROM gitpod/workspace-full

USER gitpod

# Install Conda
RUN wget -qO ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash ~/miniconda.sh -b -p ~/miniconda && \
    rm ~/miniconda.sh

# Add Conda to PATH
ENV PATH="/home/gitpod/miniconda/bin:${PATH}"
