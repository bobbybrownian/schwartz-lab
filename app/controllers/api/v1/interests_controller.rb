class Api::V1::InterestsController < ApiController

 def index
   @interests = Interest.all
   
   render json: {
     interests: @interests
   }
   
 end

 def show
  id = params[:id].to_i
   @interests = Interest.where(interest_type: id)

   render json: @interests
 end

 private


end
