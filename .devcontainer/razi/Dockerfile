#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Riccardo Vianello. All rights reserved.
# 
#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See https://go.microsoft.com/fwlink/?linkid=2090316 for license information.
#-------------------------------------------------------------------------------------------------------------

# Note: You can use any Debian/Ubuntu based image you want. 
FROM rvianello/fedora-rdkit-python:35-2021.09.3

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
    dnf install -y curl dnf-plugins-core git less libpq python3-ipython \
        python3-docutils python3-pylint python3-rstcheck \
        python3-sqlalchemy python3-sphinx python3-psycopg2 \
        twine which \
    # Install Docker CE CLI
    && dnf config-manager \
        --add-repo \
        https://download.docker.com/linux/fedora/docker-ce.repo \
    && dnf install -y docker-ce-cli \
    # Install Docker Compose
    && dnf install -y docker-compose \
    #
    # Create a non-root user to use if preferred - see https://aka.ms/vscode-remote/containers/non-root-user.
    && groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # [Optional] Add sudo support for the non-root user
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME\
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER ${USERNAME}


