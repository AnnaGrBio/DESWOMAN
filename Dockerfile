FROM ghcr.io/astral-sh/uv:bookworm-slim

RUN apt-get update && \
    apt-get upgrade --yes

# Install non-python dependencies
RUN apt-get --yes install wget ncbi-blast+ gffread tk \
    && wget http://github.com/bbuchfink/diamond/releases/download/v2.0.4/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz

# add user and set up the working directory
RUN useradd --create-home deswoman_user
USER deswoman_user
WORKDIR /home/deswoman_user

# Enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1
# Copy from the cache instead of linking since it's a mounted volume
ENV UV_LINK_MODE=copy

# Install the project's dependencies using the lockfile and settings
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev

# Installing separately from its dependencies allows optimal layer caching
COPY --chown=deswoman_user src/ pyproject.toml uv.lock README.md /home/deswoman_user/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev

# Place executables in the environment at the front of the path
ENV PATH="/home/deswoman_user/.venv/bin:$PATH"

# Reset the entrypoint, don't invoke `uv`
ENTRYPOINT []
