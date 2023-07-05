FROM gitpod/workspace-full

USER gitpod

# Install Miniconda
RUN wget -qO ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash ~/miniconda.sh -b -p ~/miniconda && \
    rm ~/miniconda.sh

# Add Conda to PATH and initialize Conda
ENV PATH="/home/gitpod/miniconda/bin:${PATH}"
RUN echo ". /home/gitpod/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
