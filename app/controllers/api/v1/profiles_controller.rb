class Api::V1::ProfilesController < ApiController

def index
  @profiles = Profile.all
  @pastmembers = Pastmember.all

  render json: {
    profiles: @profiles,
    pastmembers: @pastmembers
  }
end

end
