class Api::V1::NewstuffsController < ApiController

def index
  @newstuffs = Newstuff.all

  render json: {
    newstuffs: @newstuffs
  }
end

end
