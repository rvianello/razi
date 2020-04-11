#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Riccardo Vianello. All rights reserved.
# 
#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See https://go.microsoft.com/fwlink/?linkid=2090316 for license information.
#-------------------------------------------------------------------------------------------------------------

# Note: You can use any Debian/Ubuntu based image you want. 
FROM rvianello/fedora-rdkit-python:31-2020.03.1

# Docker Compose version
ARG COMPOSE_VERSION=1.24.0

# The image creates a non-root user with sudo access. Use the 
# "remoteUser" property in devcontainer.json to use it. On Linux, update 
# these values to ensure the container user's UID/GID matches your local values.
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Configure apt and install packages
RUN \
    # Install various tools \
    dnf install -y curl dnf-plugins-core git less postgresql which \
    # Install Docker CE CLI
    && dnf config-manager \
        --add-repo \
        https://download.docker.com/linux/fedora/docker-ce.repo \
    && dnf install -y docker-ce-cli \
    # Install Docker Compose
    # && curl -sSL "https://github.com/docker/compose/releases/download/${COMPOSE_VERSION}/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
    # && chmod +x /usr/local/bin/docker-compose
    && dnf install -y docker-compose \
    #
    # Create a non-root user to use if preferred - see https://aka.ms/vscode-remote/containers/non-root-user.
    && groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # [Optional] Add sudo support for the non-root user
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME\
    && chmod 0440 /etc/sudoers.d/$USERNAME \
    && mkdir -p /opt/venv \
    && chown -R ${USER_UID}:${USER_GID} /opt/venv

USER ${USERNAME}

RUN python3 -m venv --system-site-packages --clear /opt/venv \
    && source /opt/venv/bin/activate \
    && pip install ipython


